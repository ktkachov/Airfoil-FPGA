/*
  Mesh preparation code for the MaxCompiler implementation of Airfoil.
  The idea is to partition the mesh into chunks that will fit in the BRAM
  of the FPGA. Each of those chunks needs to be partitioned in 2, so that
  we can process one of them while reading in the other, without overlap (float buffering).
  Each of those two partitions has to be partitioned into more partitions that will be coloured
  so that any cell that is accessed will not be accessed within the next C number of accesses
  so that the arithmetic pipeline can compute the contribution of that cell to the overall value.
  The top-level chunks/partitions share halo data between them that will be reduced (with addition)
  on the host. The code has to determine the halos between every pair of partitions and schedule them
  for streaming to the chip. It also needs to compute an iteration order for the partition data. 
  The FPGA itself will not be accessing the DRAM in a random way.

  @author: Kyrylo Tkachov (kt208@imperial.ac.uk)
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "metis.h"
#include "airfoil_utils.h"

#ifdef RUN_FPGA
#include <MaxCompilerRT.h>
#define DESIGN_NAME ResCalcSim
#endif


#include "airfoil_kernels.h"

#define CELLS_PER_PARTITION (1<<13)
#define EDGES_PER_PARTITION (CELLS_PER_PARTITION * 2)
#define NODES_PER_PARTITION (CELLS_PER_PARTITION)

/*This depends on the arithmetic pipeline depth on the FPGA*/
#define PIPELINE_LATENCY 1
#define NUM_MICRO_PARTITIONS (20 * PIPELINE_LATENCY)

#define PRIME 60013
#define SMALL_PRIME 10007

#define NOP_EDGE UINT_MAX
/*colourNames is used to create a DOT file that dumps graphs in a renderable format*/

float gam;
float gm1;
float cfl;
float eps;
float qinf[4];

const char* colourNames[] = {
  "aqua",
  "olive",
  "red",
  "purple",
  "orchid",
  "orange",
  "tomato",
  "yellow",
  "cyan",
  "beige",
  "pink",
  "gray0",
  "aquamarine4"
};

/*Undirected graph structure, represented as ajacency lists*/
typedef struct graph_struct {
  uint32_t** adj_list;
  uint32_t* adj_sizes;
  uint32_t* colours;
  uint32_t num_nodes;
  uint32_t num_colours;
} ugraph;

typedef struct colour_list_struct {
  uint32_t** nodes;
  uint32_t num_colours;
  uint32_t* sizes;
} colour_list;

typedef struct size_vector_struct {
  uint32_t nodes, cells, edges, halo_cells, halo_nodes;
  uint32_t iph_cells, iph_nodes;
  uint32_t nhd1_cells, nhd1_nodes, nhd1_edges;
  uint32_t nhd2_cells, nhd2_nodes, nhd2_edges;
} size_vector_t;

typedef struct edge_address_struct {
  uint32_t node[2];
  uint32_t cell[2];
} edge_struct;

typedef struct internal_partition_struct {
  arr_t cells;
  arr_t edges;
  arr_t nodes;
  ugraph* cg; /*Connectivity graph for internal partitions*/
  uint32_t* partitionsOrdered;
  arr_t edgesOrdered;
  uint32_t nop_edges;

  edge_struct* edgeStructsOrdered;

  uint32_t* c2n; /* cells to nodes local map*/
  hash_map* g2l_nodes; /*global to local node numbers*/
  hash_map* g2l_cells; /*global to local cell numbers*/
  arr_t* parts_cells;  /* Bottom level partition cells*/
  arr_t* parts_edges; /*Bottom level partition edges globally numbered*/
  arr_t* parts_nodes; /*Bottom level partition nodes globally numbered*/
} ipartition;

typedef struct partition_struct {
  uint32_t* enodes;
  uint32_t* ecells;
  uint32_t* neighbours;
  uint32_t nneighbours;
  uint32_t* c2n; /*Cells to nodes local map*/

  ipartition iparts[3]; /*Internal partitions: 0,1 are partitions, 2 is intra partition halo, locally numbered*/

  hash_map* g2l_nodes; /*mapping of global node number to local node number*/
  hash_map* g2l_cells; /*mapping of global cell number to local cell number*/

  arr_t haloCells;
  arr_t haloNodes;
  arr_t edges;
  arr_t cells;
  arr_t non_halo_cells;
  arr_t nodes;
  arr_t* hrCells; /*Halo region cells*/
  arr_t* hrNodes;

  hash_map* nodeAddressMap;
  hash_map* cellAddressMap;
  arr_t nodesOrdered;
  arr_t haloNodesOrdered;
  arr_t cellsOrdered;
  arr_t haloCellsOrdered;
} partition;

uint32_t elem(uint32_t*, uint32_t, uint32_t);
uint32_t elemArr(arr_t*, uint32_t);
void colourGraph(ugraph*);
uint32_t removeDups(uint32_t*, uint32_t);
uint32_t min(uint32_t, uint32_t);
uint32_t max(uint32_t, uint32_t);
void showGraph(ugraph*);
/*
  Takes an array specifying the partition number of each node, a map of elements to nodes, 
  the dimension of the map (4 for cells, 2 for edges),
  the number of elements, and the number of partitions
*/
ugraph* generateGraph(uint32_t* npart, uint32_t* map, uint32_t dim, uint32_t len, uint32_t num_parts) {
  uint32_t** adj_list = malloc(num_parts * sizeof(*adj_list));
  uint32_t* adj_sizes = calloc(sizeof(*adj_sizes), num_parts);
  for (uint32_t i = 0; i < num_parts; ++i) {
    adj_list[i] = malloc(num_parts * sizeof(*adj_list[i]));
  } 
  for (uint32_t i = 0; i < len; ++i) {
    uint32_t p[dim];
    for (uint32_t j = 0; j < dim; ++j) {
      p[j] = npart[map[dim*i+j]];
    }
    for (uint32_t j = 0; j < dim; ++j) {
      for (uint32_t k = 0; k < dim; ++k) {
        if (p[j] != p[k] && !elem(adj_list[p[j]], adj_sizes[p[j]], p[k])) {
          adj_list[p[j]][adj_sizes[p[j]]++] = p[k];
        }
      }
    }
  }
  ugraph* res = malloc(sizeof(*res));
  res->adj_list = adj_list;
  res->adj_sizes = adj_sizes;
  res->num_nodes = num_parts;
  res->num_colours = 0;
  res->colours = NULL;
  return res;
}

/*
  Greedy colouring algorithm
*/
void colourGraph(ugraph* graph) {
  uint32_t nv = graph->num_nodes; //Number of vertices, also maximum colours
  uint32_t* colours = malloc(graph->num_nodes * sizeof(*colours));
  for (uint32_t i = 0; i < nv; ++i) {
    colours[i] = UINT_MAX;
  }
  uint32_t ncolours = 1;
  for (uint32_t i = 0; i < nv; ++i) {
    uint32_t c = 0;
    uint32_t repeat = 1;
    while (repeat) {
      uint32_t a = 1; //colour availble
      for (uint32_t v = 0; v < graph->adj_sizes[i] && a; ++v) {
        if (colours[graph->adj_list[i][v]] == c) {
          a = 0;
        }
      }
      if (a) {
        colours[i] = c;
        repeat = 0;
      }
      ++c;
    } 
    if (c == ncolours) {
      ++ncolours;
    }
  }
  graph->colours = colours;
  graph->num_colours = ncolours-1;
}

uint32_t** toAdjacencyMatrix(ugraph* g) {
  uint32_t** mat = malloc(g->num_nodes * sizeof(*mat));
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    mat[i] = malloc(g->num_nodes * sizeof(*mat[i]));
  }
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    for (uint32_t j = 0; j < g->num_nodes; ++j) {
      mat[i][j] = elem(g->adj_list[i], g->adj_sizes[i], j);
    }
  }
  return mat;
}

void destroyMatrix(uint32_t** mat, uint32_t s) {
  for (uint32_t i = 0; i < s; ++i) {
    free(mat[i]);
  }
  free(mat);
}

int validSchedule(ugraph* g, uint32_t* sch, uint32_t interval) {
  for (uint32_t i = 0; i < g->num_nodes-1; ++i) {
    for (uint32_t j = i+1; j < (i+1+interval) % g->num_nodes; ++j) {
      if (elem(g->adj_list[sch[i]], g->adj_sizes[sch[i]], sch[j])) {
        printf("ERROR! %d is within reach of %d\n", sch[i], sch[j]);
        return 0;
      }
    }
  }
  return 1;
}

int validPos(uint32_t* arr, uint32_t node, uint32_t nnodes, uint32_t count, uint32_t** mat, uint32_t interval) {
//  if (elem(arr, count, node)) {
//    return 0;
//  }
  if (count <= interval) {
    for (uint32_t i = 0; i < count; ++i) {
      if (mat[arr[i]][node]) {
        return 0;
      }
    }
    return 1;
  }
  if (count > nnodes - interval) {
    for (uint32_t i = 0; i < interval - (nnodes - count); ++i) {
      if (mat[arr[i]][node]) {
        return 0;
      }
    }
  }
  for (uint32_t i = count - interval - 1; i < count; ++i) {
    if (mat[arr[i]][node]) {
      return 0;
    }
  }
  return 1;
}

int sched(uint32_t n, uint32_t* res, uint32_t* count, short* elems, uint32_t** mat,  uint32_t nnodes, uint32_t interval) {
  res[(*count)++] = n;
  if ((*count) == nnodes) {
    return 1;
  }
  elems[n] = 1;
  for (uint32_t i = 0; i < nnodes; ++i) {
    if (!elems[i] && !mat[n][i] && validPos(res, i, nnodes, *count, mat, interval)) {
      if (sched(i, res, count, elems, mat, nnodes, interval)) {
        return 1;
      }
    }
  }
 // printf("returning 0 with n = %d, count=%d\n", n, *count);
  --(*count);
  elems[n] = 0;
  return 0;
}

uint32_t* scheduleGraph2(ugraph* g, uint32_t interval) {
  uint32_t* res = malloc(g->num_nodes * sizeof(*res));
  uint32_t count = 0;
  short* elems = calloc(g->num_nodes, sizeof(*elems));
  uint32_t** mat = toAdjacencyMatrix(g);
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    if (sched(i, res, &count, elems, mat,  g->num_nodes, interval)) {
      destroyMatrix(mat, g->num_nodes);
      free(elems);
      return res;
    }
  }
  printf("ERROR! could not schedule graph!\n");
  return NULL;
}

uint32_t* scheduleGraph(ugraph* g, uint32_t interval) {
  uint32_t* res = malloc(g->num_nodes * sizeof(*res));
  uint32_t** mat = toAdjacencyMatrix(g);
  short* elems = calloc(g->num_nodes, sizeof(*elems));
  stack* counts = createStack();
  uint32_t count = 0;
  stack* st = createStack();
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    stack_push(st, g->num_nodes - 1 - i);
  }
  stack_push(counts, g->num_nodes);
  while (!stack_isEmpty(st)) {
    uint32_t n = stack_pop(st);
    res[count++] = n;
    elems[n] = 1;
    if (count == g->num_nodes) {
      destroyStack(st);
      destroyStack(counts);
      destroyMatrix(mat, g->num_nodes);
      free(elems);
      return res;
    }
    uint32_t children = 0;
    for (uint32_t i = 0; i < g->num_nodes; ++i) {
      if (n != i && !mat[i][n] && !elems[i] && validPos(res, i, g->num_nodes, count, mat, interval)) {
        stack_push(st, i);
        ++children;
      }
    }
    if (children != 0) {
      stack_push(counts, children);
    } else {
        elems[n] = 0;
        --count;
        counts->top->val--;
        while (!stack_isEmpty(counts) && stack_peek(counts) == 0) {
          stack_pop(counts);
          elems[res[--count]] = 0;
          if (!stack_isEmpty(counts)) {
            counts->top->val--;
          }
        }
    }
  }
  printf("ERROR! Could not schedule graph!\n");
  free(elems);
  destroyStack(st);
  return NULL;
}

colour_list* toColourList(ugraph* g) {
  uint32_t** nodes = malloc(g->num_colours * sizeof(*nodes));
  uint32_t* sizes = calloc(g->num_colours, sizeof(*sizes));
  for (uint32_t n = 0; n < g->num_nodes; ++n) {
    ++sizes[g->colours[n]];
  }
  for (uint32_t c = 0; c < g->num_colours; ++c) {
    nodes[c] = malloc(sizes[c] * sizeof(nodes[c]));
  }
  uint32_t* counts = calloc(g->num_colours, sizeof(*counts));

  for (uint32_t n = 0; n < g->num_nodes; ++n) {
    nodes[g->colours[n]][counts[g->colours[n]]++] = n;
  }
  colour_list* res = malloc(sizeof(*res));
  res->nodes = nodes;
  res->sizes = sizes;
  res->num_colours = g->num_colours;
  free(counts);
  return res;
}

inline void printarray(float* data, uint32_t len, const char* file_name) {
  FILE* flog;
  flog = fopen(file_name, "w");
  for (uint32_t i = 0; i < len; ++i) {
    fprintf(flog, "%lf\n", data[i]);
  }
  fclose(flog);
}


/*
  Taken from the OP2 library
*/
uint32_t removeDups(uint32_t* arr, uint32_t len) {
  uint32_t i, j;
  j = 0; 
  for (i = 1; i < len; ++i) {
    if (arr[i] != arr[j]) {
      ++j;
      arr[j] = arr[i];
    }
  }
  len = (j + 1);
  return len;
}

void removeDupsArr(arr_t* a) {
  uint32_t i, j;
  j = 0;
  for (i = 1; i < a->len; ++i) {
    if (a->arr[i] != a->arr[j]) {
      ++j;
      a->arr[j] = a->arr[i];
    }
  }
  a->len = j+1;
  a->arr = realloc(a->arr, a->len * sizeof(*(a->arr)));
}

uint32_t elem(uint32_t* arr, uint32_t len, uint32_t x) {
  for (uint32_t i = 0; i < len; ++i) {
    if (arr[i] == x) {
      return 1;
    }
  }
  return 0;
}

inline uint32_t elemArr(arr_t* a, uint32_t x) {
  for (uint32_t i = 0; i < a->len; ++i) {
    if (a->arr[i] == x) {
      return 1;
    }
  }
  return 0;
}

inline uint32_t min(uint32_t a, uint32_t b) {
  return a < b ? a : b;
}

inline uint32_t max(uint32_t a, uint32_t b) {
  return a > b ? a : b;
}

void showSizeVector(size_vector_t* sv) {
  printf("{nodes:%d , cells:%d, edges:%d, haloNodes:%d, haloCells:%d, nhd1Nodes:%d, nhd1Cells:%d, ", 
        sv->nodes, sv->cells, sv->edges, sv->halo_nodes, sv->halo_cells, sv->nhd1_nodes, sv->nhd1_cells);
  printf("nhd1Edges:%d, iphNodes:%d, iphCells:%d, nhd2Nodes:%d, nhd2Cells:%d, nhd2Edges:%d}\n",
        sv->nhd1_edges, sv->iph_nodes, sv->iph_cells, sv->nhd2_nodes, sv->nhd2_cells, sv->nhd2_edges);
}

/*
  Render a partition graph into a DOT file, can specify partition structs 
  with filled in hrCells fields for edge weights or leave the 2nd argument NULL
  for general ugraph processing
*/
void generateDotGraph(ugraph* g, const char* fileName, partition* ps) {
  FILE* fp = fopen(fileName, "w");
  fprintf(fp, "graph mesh {\n");

  if (g->colours != NULL) {
    for (uint32_t i = 0; i < g->num_nodes; ++i) {
      fprintf(fp, "\t%u[color=\"%s\"];\n", i, colourNames[g->colours[i]]);
    }
  }  
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    for (uint32_t j = 0; j < g->adj_sizes[i]; ++j) {
      if (g->adj_list[i][j] > i) {
        if (ps != NULL) {
          fprintf(fp, "\t%u -- %u [label=%u];\n", i, g->adj_list[i][j], ps[i].hrCells[j].len);
        } else {
          fprintf(fp, "\t%u -- %u;\n", i, g->adj_list[i][j]);
        }
      }
    }
  }
  fprintf(fp, "}");
  fclose(fp);
}

void generateScheduleDotGraph(uint32_t* sch, uint32_t s, const char* fileName) {
  FILE* fp = fopen(fileName, "w");
  fprintf(fp, "digraph schedule{\n");
  fprintf(fp, "layout=neato;\n");
  fprintf(fp, "%d[color=\"red\"];\n", sch[0]);
  fprintf(fp, "%d -> %d;\n", sch[s-1], sch[0]);
  for (uint32_t i = 0; i < s-1; ++i) {
    fprintf(fp, "%d -> %d;\n", sch[i], sch[i+1]);
  }
  fprintf(fp, "}");
  fclose(fp);
}

void showColourListSizes(colour_list* cl) {
  for (uint32_t c = 0; c < cl->num_colours; ++c) {
    printf("colour %d has %d elements\n", c, cl->sizes[c]);
  }
}


/*Dump an undirected graph to stdout*/
void showGraph(ugraph* g) {
  for (uint32_t i = 0; i < g->num_nodes; ++i) {
    if (g->colours != NULL) {
      printf("%u (colour %d): ", i, g->colours[i]);
    } else {
      printf("%u :", i);
    }
    for (uint32_t j = 0 ; j < g->adj_sizes[i]; ++j) {
      printf("%u, ", g->adj_list[i][j]);
    }
    printf("\n");
  }
  if (g->colours != NULL) {
    printf("number of colours: %d\n", g->num_colours);
  }
}

int main(int argc, char* argv[]) {

  uint32_t *becell, *ecell, *bound, *bedge, *edge, *cell;
  float *x, *q, *qold, *adt, *res;

  uint32_t nnode,ncell,nedge,nbedge;

 
  /* read in grid */

  printf("reading in grid \n");
  FILE *fp;
  if ( (fp = fopen("./new_grid.dat","r")) == NULL) { ///new_grid.dat
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%u %u %u %u \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell = malloc(4*ncell*sizeof(*cell));
  edge = malloc(2*nedge*sizeof(*edge));
  ecell = malloc(2*nedge*sizeof(*ecell));
  bedge =  malloc(2*nbedge*sizeof(*bedge));
  becell = malloc( nbedge*sizeof(*becell));
  bound =  malloc( nbedge*sizeof(*bound));

  x =  malloc(2*nnode*sizeof(*x));
  q =  malloc(4*ncell*sizeof(*q));
  qold = malloc(4*ncell*sizeof(*qold));
  res = malloc(4*ncell*sizeof(*res));
  adt = malloc( ncell*sizeof(*adt));

  for (uint32_t n=0; n<nnode; n++) {
    if (fscanf(fp,"%f %f \n",&x[2*n], &x[2*n+1]) != 2) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (uint32_t n=0; n<ncell; n++) {
    if (fscanf(fp,"%u %u %u %u \n",&cell[4*n ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (uint32_t n=0; n<nedge; n++) {
    if (fscanf(fp,"%u %u %u %u \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (uint32_t n=0; n<nbedge; n++) {
    if (fscanf(fp,"%u %u %u %u \n",&bedge[2*n],&bedge[2*n+1],
                                   &becell[n], &bound[n]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  fclose(fp);

  clock_t start, end;
  start = clock();
 
  /* Cell pointer for METIS. Since all cells have 4 nodes, this is not very interesting, but METIS requires it*/
  uint32_t* cptr = malloc((ncell+1) * sizeof(*cptr));
  for (uint32_t i = 0; i < ncell+1; ++i) {
    cptr[i] =4*i;
  }
  int objval;
  uint32_t *cpart = malloc(ncell * sizeof(*cpart));
  uint32_t *npart = malloc(nnode * sizeof(*npart));
  int num_parts = ncell / CELLS_PER_PARTITION;
  int nn = nnode;
  int nc = ncell;
  printf("Partitioning %u cells in %u partitions...\n", ncell, num_parts);
  /*Initial partitioning of the mesh by cells*/
  METIS_PartMeshNodal(&nc, &nn, (int*)cptr, (int*)cell, NULL, NULL, &num_parts, NULL, NULL, &objval, (int*)cpart, (int*)npart);

  /*Initialising our partition structures*/
  partition* ps = malloc(num_parts * sizeof(*ps));
  
  hash_set* edge_maps[num_parts];
  hash_set* node_maps[num_parts];
  hash_set* cell_maps[num_parts];
  hash_set* halo_cells[num_parts];
  for (uint32_t i = 0; i < num_parts; ++i) {
    node_maps[i] = createHashSet(PRIME);
    edge_maps[i] = createHashSet(PRIME);
    cell_maps[i] = createHashSet(PRIME);
    halo_cells[i] = createHashSet(SMALL_PRIME);
    initArr(&ps[i].edges, EDGES_PER_PARTITION);
    initArr(&ps[i].nodes, NODES_PER_PARTITION);
    initArr(&ps[i].cells, CELLS_PER_PARTITION);
    initArr(&ps[i].haloNodes, NODES_PER_PARTITION / 10);
    initArr(&ps[i].haloCells, CELLS_PER_PARTITION / 10);
    initArr(&ps[i].non_halo_cells, CELLS_PER_PARTITION);
  }

  for (uint32_t i = 0; i < nedge; ++i) {
    uint32_t n1 = npart[edge[2*i]];
//    uint32_t n2 = npart[edge[2*i+1]];
    addToHashSet(edge_maps[n1], i);

    uint32_t c1 = ecell[2*i];
    uint32_t c2 = ecell[2*i+1];

    addToHashSet(cell_maps[n1], c1);
    addToHashSet(cell_maps[n1], c2);

    for (short j = 0; j < 4; ++j) {
      if (npart[cell[4*c1+j]] != n1) {
        addToHashSet(halo_cells[n1], c1);
      }
      if (npart[cell[4*c2+j]] != n1) {
        addToHashSet(halo_cells[n1], c2);
      }
      addToHashSet(node_maps[n1], cell[4*c1 + j]);
      addToHashSet(node_maps[n1], cell[4*c2 + j]);
    }
  }
  for (uint32_t i = 0; i < num_parts; ++i) {
    ps[i].edges = *toArr(edge_maps[i]);
    ps[i].cells = *toArr(cell_maps[i]);
    ps[i].nodes = *toArr(node_maps[i]);
    ps[i].haloCells = *toArr(halo_cells[i]);

    hash_set* diff = setDiff(cell_maps[i], halo_cells[i]);
    ps[i].non_halo_cells = *toArr(diff);
    destroyHashSet(diff);
  //  destroyHashSet(halo_cells[i]);
  }
  printf("Assigned edges, cells, nodes and halo cells to partitions\n");
  
  /*Clean up the structures populated above and initialise the local cells-to-nodes maps*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    ps[i].c2n = malloc(4 * ps[i].cells.len * sizeof(*ps[i].c2n));
  }


  uint32_t num_nop_edges = 0;
  uint32_t p_edges = 0;
  for (uint32_t i = 0; i < num_parts; ++i) {
    printf("Partition %d has %d cells, %d nodes, %d edges\n", i, ps[i].cells.len, ps[i].nodes.len, ps[i].edges.len);
    /*
      We use hash maps to store the mapping of global cell and node numbers to local numbers.
      We need local cell and node numbers in order to partition the partitions internally.
    */
    ps[i].g2l_nodes = createHashMap(PRIME);
    ps[i].g2l_cells = createHashMap(PRIME);

    for (uint32_t j = 0; j < ps[i].nodes.len; ++j) {
      addToHashMap(ps[i].g2l_nodes, ps[i].nodes.arr[j], j);
    }
    for (uint32_t j = 0; j < ps[i].cells.len; ++j) {
      addToHashMap(ps[i].g2l_cells, ps[i].cells.arr[j], j);
    }
    for (uint32_t j = 0; j < ps[i].cells.len; ++j) {
      for (short k = 0; k < 4; ++k) {
        uint32_t v = getValue(ps[i].g2l_nodes, cell[4*ps[i].cells.arr[j] + k]);
        if (v == EMPTY_ENTRY) {
          printf("Error! ps[%d].c2n[%d] is a node not in node set!\n", i, v);
        }
        ps[i].c2n[4*j+k] = v;
      }
    }
    /*Preparations for partitioning the partition into 2 partitions*/
    uint32_t* pcptr = malloc((ps[i].cells.len + 1) * sizeof(*pcptr));
    for (uint32_t j = 0; j < ps[i].cells.len + 1; ++j) {
      pcptr[j] = 4*j;
    }
    uint32_t* pcpart = malloc(ps[i].cells.len * sizeof(*pcpart));
    uint32_t* pnpart = malloc(ps[i].nodes.len * sizeof(*pnpart));
    uint32_t nparts = 2;
    METIS_PartMeshNodal((int*)&ps[i].cells.len, (int*)&ps[i].nodes.len, (int*)pcptr, (int*)ps[i].c2n, NULL, NULL,(int*) &nparts, NULL, NULL, &objval, (int*)pcpart, (int*)pnpart);
    free(pcptr);
    pcptr = NULL;

    /*Here we determine the two internal partitions and the intra-partition halo cells and edges*/
    for (short j = 0; j < 3; ++j) {
      initArr(&ps[i].iparts[j].cells, CELLS_PER_PARTITION / 2);
      initArr(&ps[i].iparts[j].edges, EDGES_PER_PARTITION / 2);
      initArr(&ps[i].iparts[j].nodes, NODES_PER_PARTITION / 2);
    }
    hash_set* edge_sets[3];
    hash_set* cell_sets[3];
    hash_set* node_sets[3];
    for (short j = 0; j < 3; ++j) {
      edge_sets[j] = createHashSet(SMALL_PRIME);
      cell_sets[j] = createHashSet(SMALL_PRIME);
      node_sets[j] = createHashSet(SMALL_PRIME);
    }
    for (uint32_t j = 0; j < ps[i].edges.len; ++j) {
      
      uint32_t n1 = pnpart[getValue(ps[i].g2l_nodes, edge[2*ps[i].edges.arr[j]])];

      addToHashSet(edge_sets[n1], ps[i].edges.arr[j]);

      uint32_t c1 = getValue(ps[i].g2l_cells, ecell[2*ps[i].edges.arr[j]]);
      uint32_t c2 = getValue(ps[i].g2l_cells, ecell[2*ps[i].edges.arr[j] + 1]);
      if (c1 == EMPTY_ENTRY || c2 == EMPTY_ENTRY) {
        printf("ERROR! cell hash map does not contain cell values!\n");
      }
      addToHashSet(cell_sets[n1], c1);
      addToHashSet(cell_sets[n1], c2);

      for (short k = 0; k < 4; ++k) {
        uint32_t node1 = getValue(ps[i].g2l_nodes, cell[4*ps[i].cells.arr[c1] + k]);
        uint32_t node2 = getValue(ps[i].g2l_nodes, cell[4*ps[i].cells.arr[c2] + k]);

        if (node1 == EMPTY_ENTRY || node2 == EMPTY_ENTRY) {
          printf("ERROR! node hash map does not contain node values!\n");
        }

        if (pnpart[node1] != n1) {
          addToHashSet(cell_sets[2], c1);
        }
        if (pnpart[node2] != n1) {
          addToHashSet(cell_sets[2], c2);
        }
        addToHashSet(node_sets[n1], node1);
        addToHashSet(node_sets[n1], node2);
      }

    }
    destroyHashSet(node_sets[2]);
    node_sets[2] = setIntersection(node_sets[0], node_sets[1]);
    for (short j = 0; j < 3; ++j) {
      ps[i].iparts[j].edges = *toArr(edge_sets[j]);
      ps[i].iparts[j].cells = *toArr(cell_sets[j]);
      ps[i].iparts[j].nodes = *toArr(node_sets[j]);
    }

    for (uint32_t j = 0; j < 2; ++j) {
      ipartition* ip = &ps[i].iparts[j];
      ip->g2l_nodes = createHashMap(SMALL_PRIME);
      ip->g2l_cells = createHashMap(SMALL_PRIME);
      ip->c2n = malloc(ip->cells.len * 4 * sizeof(*ip->c2n));

      //printf("Subpartition %d of partition %d has %d nodes, %d edges, %d cells\n", j, i, ip->nodes.len, ip->edges.len, ip->cells.len);

      for (uint32_t k = 0; k < ip->nodes.len; ++k) {
        addToHashMap(ip->g2l_nodes, ip->nodes.arr[k], k);
      }
      for (uint32_t k = 0; k < ip->cells.len; ++k) {
        addToHashMap(ip->g2l_cells, ip->cells.arr[k], k);
      }

      for (uint32_t k = 0; k < ip->cells.len; ++k) {
        for (short kk = 0; kk < 4; ++kk) {
          uint32_t n = getValue(ip->g2l_nodes, ps[i].c2n[4*ip->cells.arr[k] + kk]);
          ip->c2n[4*k+kk] = n;
          if (n == EMPTY_ENTRY) {
            printf("Error! map @ %d does not contain node\n", k);
          }
        }
      }

      pcptr = malloc((ip->cells.len + 1) * sizeof(*pcptr));
      for (uint32_t k = 0; k < ip->cells.len + 1; ++k) {
        pcptr[k] = 4*k;
      }
      nparts = NUM_MICRO_PARTITIONS;
      uint32_t* intcpart = malloc((ip->cells.len+1) * sizeof(*intcpart));
      uint32_t* intnpart = malloc(ip->nodes.len * sizeof(*intnpart));
      METIS_PartMeshNodal((int*)&ip->cells.len,(int*)&ip->nodes.len, (int*)pcptr, (int*)ip->c2n, NULL, NULL, (int*)&nparts, NULL, NULL, &objval, (int*)intcpart, (int*)intnpart);
      free(pcptr);

      ip->cg = generateGraph(intnpart, ip->c2n, 4, ip->cells.len, nparts);
      colourGraph(ip->cg);
      ip->parts_nodes = malloc(nparts * sizeof(*ip->parts_nodes));
      ip->parts_cells = malloc(nparts * sizeof(*ip->parts_cells));
      ip->parts_edges = malloc(nparts * sizeof(*ip->parts_edges));
      hash_set* bottom_edge_sets[nparts];
      for (uint32_t k = 0; k < nparts; ++k) {
        bottom_edge_sets[k] = createHashSet(SMALL_PRIME);
      }
      for (uint32_t k = 0; k < ip->edges.len; ++k) {
        uint32_t n[2];
        n[0] = getValue(ip->g2l_nodes, getValue(ps[i].g2l_nodes, edge[2*ip->edges.arr[k]]));
        n[1] = getValue(ip->g2l_nodes, getValue(ps[i].g2l_nodes, edge[2*ip->edges.arr[k] + 1]));
        addToHashSet(bottom_edge_sets[intnpart[n[0]]], ip->edges.arr[k]);
      }
      uint32_t se = 0;
      for (uint32_t k = 0; k < nparts; ++k) {
        ip->parts_edges[k] = *toArr(bottom_edge_sets[k]);
        destroyHashSet(bottom_edge_sets[k]);
        se += ip->parts_edges[k].len;
        //printf("bottom level partition %d of partition %d of partition %d has %d edges\n", k, j, i, ip->parts_edges[k].len);
      }
      p_edges += se;
 //     printf("internal graph for partition %d of partition %d is:\n", j, i);
 //     showGraph(ip->cg);
    }
 
    for (short k = 0; k < 2; ++k) {
      hash_set* newCells = setDiff(cell_sets[k], cell_sets[2]);
      destroyHashSet(cell_sets[k]);
      ps[i].iparts[k].cells = *toArr(newCells); 
      destroyHashSet(newCells);

      hash_set* newNodes = setDiff(node_sets[k], node_sets[2]);
      destroyHashSet(node_sets[k]);
      ps[i].iparts[k].nodes = *toArr(newNodes);
      destroyHashSet(newNodes);
    }
/*
    for (short k = 0; k < 3; ++k) {
      printf("sub-partition %d of partition %d after common elements removal has %d cells, %d nodes, %d edges\n", 
              k,
              i,
              ps[i].iparts[k].cells.len,
              ps[i].iparts[k].nodes.len,
              ps[i].iparts[k].edges.len
            );
    }
*/

  }

  uint32_t part_nodes = 0;
  uint32_t part_cells = 0;
  uint32_t part_edges = 0;
  for (uint32_t i = 0; i < num_parts; ++i) {
    part_nodes += ps[i].nodes.len;
    part_cells += ps[i].cells.len;
    part_edges += ps[i].edges.len;
    printf("partition %u has: %u edges, %u nodes, %u cells, %u halo cells\n", i, ps[i].edges.len, ps[i].nodes.len, ps[i].non_halo_cells.len, ps[i].haloCells.len);
  }
  printf("Counted %d nodes, %d cells, %d edges\n", part_nodes, part_cells, part_edges);
  printf("Generating top level partition graph\n");
  ugraph* pg = generateGraph(npart, cell, 4, ncell, num_parts);
  if (!pg) {
    printf("ERROR! partition graph is NULL!\n");
    return 1;
  }


  /*Initialise the per-neighbour halo regions for each partition*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    ps[i].neighbours = malloc(pg->adj_sizes[i] * sizeof(*(ps[i].neighbours)));
    ps[i].nneighbours = pg->adj_sizes[i];
    memcpy(ps[i].neighbours, pg->adj_list[i], pg->adj_sizes[i] * sizeof(*(ps[i].neighbours)));
    ps[i].hrCells = malloc(ps[i].nneighbours * sizeof(*ps[i].hrCells));
    ps[i].hrNodes = malloc(ps[i].nneighbours * sizeof(*ps[i].hrNodes));
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      initArr(&ps[i].hrCells[j], 64);
    }
  }

  printf("Initialised halo regions\n");

  /*Determine the cells that are shared between every pair of partitions*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    hash_set* hr_cell_sets[ps[i].nneighbours];
    hash_set* halo_nodes = createHashSet(101);
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      hr_cell_sets[j] = setIntersection(halo_cells[i], halo_cells[ps[i].neighbours[j]]);
      ps[i].hrCells[j] = *toArr(hr_cell_sets[j]);
      destroyHashSet(hr_cell_sets[j]);

      hash_set* nodesIntersection = setIntersection(node_maps[i], node_maps[ps[i].neighbours[j]]);
      ps[i].hrNodes[j] = *toArr(nodesIntersection);
      hash_set* nodes_temp = setUnion(halo_nodes, nodesIntersection);
      destroyHashSet(halo_nodes);
      halo_nodes = nodes_temp;
      destroyHashSet(nodesIntersection);
    }
    ps[i].haloNodes = *toArr(halo_nodes);
    destroyHashSet(halo_nodes);
  }

  for (uint32_t i = 0; i < num_parts; ++i) {
    initArr(&ps[i].nodesOrdered, NODES_PER_PARTITION);
    initArr(&ps[i].haloNodesOrdered, NODES_PER_PARTITION / 100);
    initArr(&ps[i].cellsOrdered, CELLS_PER_PARTITION);
    initArr(&ps[i].haloCellsOrdered, CELLS_PER_PARTITION / 100);
    ps[i].nodeAddressMap = createHashMap(SMALL_PRIME);
    ps[i].cellAddressMap = createHashMap(SMALL_PRIME);
    short a[3] = {0, 2, 1};
    for (short p = 0; p < 3; ++p) {
      for (uint32_t n = 0; n < ps[i].iparts[a[p]].nodes.len; ++n) {
        uint32_t node = ps[i].nodes.arr[ps[i].iparts[a[p]].nodes.arr[n]];
        if (!elemArr(&ps[i].haloNodes, node)) {
          addToArr(&ps[i].nodesOrdered, node);
        } else {
          addToArr(&ps[i].haloNodesOrdered, node);
        }
      }

      for (uint32_t c = 0; c < ps[i].iparts[a[p]].cells.len; ++c) {
        uint32_t ctemp = ps[i].cells.arr[ps[i].iparts[a[p]].cells.arr[c]];
        if (!elemArr(&ps[i].haloCells, ctemp)) {
          addToArr(&ps[i].cellsOrdered, ctemp);
        } else {
          addToArr(&ps[i].haloCellsOrdered, ctemp);
        }
      }
    }
    uint32_t node_count = 0;
    for (uint32_t n = 0; n < ps[i].nodesOrdered.len; ++n) {
      node_count += addToHashMap(ps[i].nodeAddressMap, ps[i].nodesOrdered.arr[n], node_count);
    }
    for (uint32_t n = 0; n < ps[i].haloNodesOrdered.len; ++n) {
      node_count += addToHashMap(ps[i].nodeAddressMap, ps[i].haloNodesOrdered.arr[n], node_count);
    }
    uint32_t cell_count = 0;
    for (uint32_t n = 0; n < ps[i].cellsOrdered.len; ++n) {
      cell_count += addToHashMap(ps[i].cellAddressMap, ps[i].cellsOrdered.arr[n], cell_count);
    }
    for (uint32_t n = 0; n < ps[i].haloCellsOrdered.len; ++n) {
      cell_count += addToHashMap(ps[i].cellAddressMap, ps[i].haloCellsOrdered.arr[n], cell_count);
    }
    //printf("Partition %d has %d ordered nodes and %d ordered halo nodes, total nodes in partition: %d\n", i, ps[i].nodesOrdered.len, ps[i].haloNodesOrdered.len, ps[i].nodes.len);

    for (short p = 0; p < 2; ++p) { 
      initArr(&ps[i].iparts[p].edgesOrdered, EDGES_PER_PARTITION / 2);
      ugraph* g = ps[i].iparts[p].cg;
    //  colour_list* cl = toColourList(g);
      ps[i].iparts[p].partitionsOrdered = scheduleGraph2(g, PIPELINE_LATENCY);
      printf("Scheduled partition graph for partition %d of %d is: \n", p, i);
      printf("is schedule valid? %d\n", validSchedule(g, ps[i].iparts[p].partitionsOrdered , PIPELINE_LATENCY));
      for (uint32_t j = 0; j < g->num_nodes; ++j) {
        printf("%d, ", ps[i].iparts[p].partitionsOrdered [j]);
      }
      printf("\n");
    //  showGraph(g);


      uint32_t* partsOrdered = ps[i].iparts[p].partitionsOrdered;
      uint32_t hMax = 0;
      for (uint32_t j = 0; j < NUM_MICRO_PARTITIONS; ++j) {
        hMax = ps[i].iparts[p].parts_edges[j].len > hMax ? ps[i].iparts[p].parts_edges[j].len : hMax;
      }
      ps[i].iparts[p].nop_edges = 0;
      for (uint32_t h = 0; h < hMax; ++h) {
        for (uint32_t pp = 0; pp < NUM_MICRO_PARTITIONS; ++pp) {
          if (ps[i].iparts[p].parts_edges[partsOrdered[pp]].len <= h) {
            addToArr(&ps[i].iparts[p].edgesOrdered, NOP_EDGE);
            ps[i].iparts[p].nop_edges++;
          } else {
            addToArr(&ps[i].iparts[p].edgesOrdered, ps[i].iparts[p].parts_edges[partsOrdered[pp]].arr[h]);
          }
        }
      }
      num_nop_edges += ps[i].iparts[p].nop_edges;
      printf("Added %d nop edges in part %d of part %d, nop to edge ratio: %.4f\n", ps[i].iparts[p].nop_edges, p, i, (float)ps[i].iparts[p].nop_edges / ps[i].iparts[p].edges.len);
     
      ps[i].iparts[p].edgeStructsOrdered = malloc(ps[i].iparts[p].edgesOrdered.len * sizeof(*ps[i].iparts[p].edgeStructsOrdered));
      for (uint32_t e = 0; e < ps[i].iparts[p].edgesOrdered.len; ++e) {
        if (ps[i].iparts[p].edgesOrdered.arr[e] == NOP_EDGE) {
          ps[i].iparts[p].edgeStructsOrdered[e].node[0] = 0;
          ps[i].iparts[p].edgeStructsOrdered[e].node[1] = 0;
          ps[i].iparts[p].edgeStructsOrdered[e].cell[0] = 0;
          ps[i].iparts[p].edgeStructsOrdered[e].cell[1] = 0;
        } else {
          uint32_t ed = ps[i].iparts[p].edgesOrdered.arr[e];
          uint32_t node_g[2]; /*Globally numbered nodes*/
          uint32_t cell_g[2];  /*Globally numbered cells*/
          node_g[0] = edge[2*ed];
          node_g[1] = edge[2*ed+1];
          cell_g[0] = ecell[2*ed];
          cell_g[1] = ecell[2*ed+1];
          for (short ii = 0; ii < 2; ++ii) {
            ps[i].iparts[p].edgeStructsOrdered[e].node[ii] = getValue(ps[i].nodeAddressMap, node_g[ii]);
            ps[i].iparts[p].edgeStructsOrdered[e].cell[ii] = getValue(ps[i].cellAddressMap, cell_g[ii]);
          }
        }
        /*
        edge_struct* estr = &ps[i].iparts[p].edgeStructsOrdered[e];
        printf("Edge %d is: {%d, %d, %d, %d}\n", e, estr->node[0], estr->node[1], estr->cell[0], estr->cell[1]);
        */
      }
    printf("\n-------------------\n");
    }

 }

  uint32_t* glPartsSched = scheduleGraph(pg, 1);


  arr_t globalNodesScheduled;
  arr_t globalCellsScheduled;
  arr_t globalHaloNodesScheduled;
  arr_t globalHaloCellsScheduled;
  initArr(&globalNodesScheduled, nnode);
  initArr(&globalCellsScheduled, ncell);
  initArr(&globalHaloNodesScheduled, num_parts * NODES_PER_PARTITION / 20);
  initArr(&globalHaloCellsScheduled, num_parts * CELLS_PER_PARTITION / 20);

  size_vector_t* size_vectors = malloc(num_parts * sizeof(*size_vectors));

  edge_struct* globalEdgeStructsScheduled = malloc((nedge + num_nop_edges) * sizeof(*globalEdgeStructsScheduled));
  uint32_t* globalEdgesScheduled = malloc((nedge + num_nop_edges) * sizeof(*globalEdgesScheduled));
  uint32_t schEdges = 0;
  for (uint32_t i = 0; i < num_parts; ++i) {
    uint32_t p = glPartsSched[i];
    size_vectors[p].nodes = ps[p].nodes.len;
    size_vectors[p].cells = ps[p].cells.len;
    size_vectors[p].edges = ps[p].iparts[0].edgesOrdered.len + ps[p].iparts[1].edgesOrdered.len;
    size_vectors[p].halo_nodes = ps[p].haloNodes.len;
    size_vectors[p].halo_cells = ps[p].haloCells.len;
    size_vectors[p].iph_cells = ps[p].iparts[2].cells.len;
    size_vectors[p].iph_nodes = ps[p].iparts[2].nodes.len;
    size_vectors[p].nhd1_nodes = ps[p].iparts[0].nodes.len;
    size_vectors[p].nhd1_nodes = ps[p].iparts[0].nodes.len;
    size_vectors[p].nhd1_cells = ps[p].iparts[0].cells.len;
    size_vectors[p].nhd1_edges = ps[p].iparts[0].edgesOrdered.len;
    size_vectors[p].nhd2_nodes = ps[p].iparts[1].nodes.len;
    size_vectors[p].nhd2_cells = ps[p].iparts[1].cells.len;
    size_vectors[p].nhd2_edges = ps[p].iparts[1].edgesOrdered.len;
    /*showSizeVector(&size_vectors[p]);*/
    for (uint32_t n = 0; n < ps[p].nodesOrdered.len; ++n) {
      addToArr(&globalNodesScheduled, ps[p].nodesOrdered.arr[n]);
    }
    for (uint32_t c = 0; c < ps[p].cellsOrdered.len; ++c) {
      addToArr(&globalCellsScheduled, ps[p].cellsOrdered.arr[c]);
    }
    for (short ip = 0; ip < 2; ++ip) {
      for (uint32_t e = 0; e < ps[p].iparts[ip].edgesOrdered.len; ++e) {
        globalEdgeStructsScheduled[schEdges] = ps[p].iparts[ip].edgeStructsOrdered[e];
        globalEdgesScheduled[schEdges] = ps[p].iparts[ip].edgesOrdered.arr[e];
        ++schEdges;
      }
    }
    for (uint32_t hn = 0; hn < ps[p].haloNodesOrdered.len; ++hn) {
     // printf("adding %d at %d to globalHaloNodesScheduled\n", ps[p].haloNodesOrdered.arr[hn], globalHaloNodesScheduled.len);
     // printf("the maxsize of halonodes ordered is:%d\n", globalHaloNodesScheduled.maxLen);
      addToArr(&globalHaloNodesScheduled, ps[p].haloNodesOrdered.arr[hn]);
    }
    for (uint32_t hc = 0; hc < ps[p].haloCellsOrdered.len; ++hc) {
      addToArr(&globalHaloCellsScheduled, ps[p].haloCellsOrdered.arr[hc]);
    }
  }

   /*Diagnostic messages, etc...*/
  printf("Calculating halo regions\n");
  uint32_t total_halo_cells = 0;
  for (uint32_t i = 0; i < num_parts; ++i) {
    uint32_t total = 0;
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      printf("Halo region %u-%u has %u cells and %d nodes\n", i, ps[i].neighbours[j], ps[i].hrCells[j].len, ps[i].hrNodes[j].len);
      total += ps[i].hrCells[j].len;
    }
    total_halo_cells += total;
    printf("Total: %u halo cells and %d halo nodes\n", total, ps[i].haloNodes.len);
    printf("----------------------------------\n");
  }






  float h2nhCells = (float)total_halo_cells / (ncell - total_halo_cells);
  printf("Colouring partition graph...\n");
  colourGraph(pg);
  printf("Top level graph of mesh partitions:\n");
  showGraph(pg);
  printf("The top-level scheduled partition graph is:\n");
  for (uint32_t i = 0; i < pg->num_nodes; ++i) {
    printf("%d, ", glPartsSched[i]);
  }
  printf("\n");

  //colour_list* cl_global = toColourList(pg);
  //printf("Colour list for global partitioning graph:\n");
  const char* fileName = "meshColoured.dot";
  printf("Writing partition graph to %s ...\n", fileName);
  generateDotGraph(pg, fileName, ps);
 
   const char* schFileName = "meshSchedule.dot";
  printf("Writing schedule graph to %s ...\n", schFileName);
  generateScheduleDotGraph(glPartsSched, num_parts, schFileName);

  end = clock();
  printf("Nodes: %u, Edges: %u, Cells: %u, number of partitions: %u\n", nnode, nedge, ncell, num_parts);
  printf("ratio of halo cells to non-halo cells is: %.2f, total halo cells: %d\n", h2nhCells, total_halo_cells);
  printf("Global scheduled nodes:%d, cells: %d\n", globalNodesScheduled.len, globalCellsScheduled.len);
  printf("Global scheduled halo nodes:%d, cells:%d\n", globalHaloNodesScheduled.len, globalHaloCellsScheduled.len);
  printf("Created global schedules, total edges = %d\n", schEdges);
  printf("ratio of nop edges to edges %.2f, with %d nop edges\n", (float)num_nop_edges / nedge, num_nop_edges);
  printf("time taken: %f seconds\n", (float)(end - start)/CLOCKS_PER_SEC);


  
  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;


  float mach = 0.4f;
//   float alpha = 3.0f*atan(1.0f)/45.0f;
  float p = 1.0f;
  float r = 1.0f;
  float u = sqrtf(gam*p/r)*mach;
  float e = p/(r*gm1) + 0.5f*u*u;

  qinf[0] = r;
  qinf[1] = r*u;
  qinf[2] = 0.0f;
  qinf[3] = r*e;

  for (long n=0; n<ncell; n++) {
    for (long m=0; m<4; m++) {
        q[4*n+m] = qinf[m];
      res[4*n+m] = 0.0f;
      
    }
  }

  printf("Running save_soln and adt_calc on host to prepare data arrays...\n");  
  for (uint32_t i = 0; i < ncell; ++i) {
    save_soln(&q[4*i], &qold[4*i]);
  }
  
  for (uint32_t i = 0; i < ncell; ++i) {
    adt_calc(&x[2*cell[4*i]], &x[2*cell[4*i+1]], &x[2*cell[4*i+2]], &x[2*cell[4*i+3]], &q[4*i], &adt[i]);
  }
  
  printf("Scheduling data arrays...\n\n");
  float* q_scheduled = malloc(globalCellsScheduled.len * 4 * sizeof(*q_scheduled));
  float* adt_scheduled = malloc(globalCellsScheduled.len * sizeof(*adt_scheduled));
  for (uint32_t i = 0; i < globalCellsScheduled.len; ++i) {
    for (short j = 0; j < 4; ++j) {
      q_scheduled[4*i+j] = q[4*globalCellsScheduled.arr[i] + j];
    }
    adt_scheduled[i] = adt[globalCellsScheduled.arr[i]];
  }
  float* x_scheduled = malloc(globalNodesScheduled.len * 2 * sizeof(*x_scheduled));
  for (uint32_t i = 0; i < globalNodesScheduled.len; ++i) {
    x_scheduled[2*i] = x[2*globalNodesScheduled.arr[i]];
    x_scheduled[2*i+1] = x[2*globalNodesScheduled.arr[i]+1];
  }

  float* halo_q_scheduled = malloc(globalHaloCellsScheduled.len * 4 * sizeof(*halo_q_scheduled));
  float* halo_adt_scheduled = malloc(globalHaloCellsScheduled.len * sizeof(*halo_q_scheduled));
  for (uint32_t i = 0; i < globalHaloCellsScheduled.len; ++i) {
    for (short j = 0; j < 4; ++j) {
      halo_q_scheduled[4*i+j] = q[4*globalHaloCellsScheduled.arr[i] + j];
    }
   halo_adt_scheduled[i] = adt[globalHaloCellsScheduled.arr[i]];
  }
  float* halo_x_scheduled = malloc(globalHaloNodesScheduled.len * 2 * sizeof(*halo_x_scheduled));
  for (uint32_t i = 0; i < globalHaloNodesScheduled.len; ++i) {
    halo_x_scheduled[2*i] = x[2*globalHaloNodesScheduled.arr[i]];
    halo_x_scheduled[2*i+1] = x[2*globalHaloNodesScheduled.arr[i]+1];
  }

  #ifdef RUN_FPGA
  short isSimulation = 1;
  char* device_name = isSimulation ? "sim:sim0" : "/dev/maxeler0";
  max_maxfile_t* maxfile;
  max_device_handle_t* device;

  printf("Initialising maxfile...\n");
  maxfile = max_maxfile_init_ResCalc();
  printf("Opening device %s ... \n", device_name);
  device = max_open_device(maxfile, device_name);
  max_set_terminate_on_error(device);


  printf("Closing device %s ... \n", device_name);
  max_close_device(device);
  printf("Destroying maxfile ... \n");
  max_destroy(maxfile);
  #endif 


  
}
