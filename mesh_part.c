/*
  Mesh preparation code for the MaxCompiler implementation of Airfoil.
  The idea is to partition the mesh into chunks that will fit in the BRAM
  of the FPGA. Each of those chunks needs to be partitioned in 2, so that
  we can process one of them while reading in the other, without overlap (double buffering).
  Each of those two partitions has to be partitioned into more partitions that will be coloured
  so that any cell that is accessed will not be accessed within the next C number of accesses
  so that the arithmetic pipeline can compute the contribution of that cell to the overall value.
  The top-level chunks/partitions share halo data between them that will be reduced (with addition)
  on the host. The code has to determine the halos between every pair of partitions and schedule them
  for streaming to the chip. It also needs to compute an iteration order for the partition data. The
  FPGA itself will not be accessing the DRAM in a random way.

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

#define CELLS_PER_PARTITION (1<<14)
#define EDGES_PER_PARTITION (CELLS_PER_PARTITION * 2)
#define NODES_PER_PARTITION (CELLS_PER_PARTITION)

/*This depends on the arithmetic pipeline depth on the FPGA*/
#define BOTTOM_LEVEL_PARTITIONS 17

#define PRIME 60013
#define SMALL_PRIME 10007

/*colourNames is used to create a DOT file that dumps graphs in a renderable format*/
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
} ugraph;


typedef struct internal_partition_struct {
  arr_t cells;
  arr_t edges;
  ugraph* cg; /*Connectivity graph for internal partitions*/
  uint32_t* c2n; /* cells to nodes local map*/
  hash_map* g2l_nodes; /*global to local node numbers*/
  hash_map* g2l_cells; /*global to local cell numbers*/
  hash_map* l2g_nodes;
  hash_map* l2g_cells;
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
  hash_map* l2g_nodes;
  hash_map* l2g_cells;

  arr_t haloCells;
  arr_t haloNodes;
  arr_t edges;
  arr_t cells;
  arr_t non_halo_cells;
  arr_t nodes;
  arr_t* hrCells; /*Halo region cells*/
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
  ugraph* res = (ugraph*)malloc(sizeof(ugraph));
  res->adj_list = adj_list;
  res->adj_sizes = adj_sizes;
  res->num_nodes = num_parts;
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
}

inline void printarray(double* data, uint32_t len, const char* file_name) {
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

inline uint32_t elem(uint32_t* arr, uint32_t len, uint32_t x) {
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
}

int main(int argc, char* argv[]) {

  uint32_t *becell, *ecell, *bound, *bedge, *edge, *cell;
  double *x, *q, *qold, *adt, *res;

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

  /*
    Cellse is the cells-to-edges map, cellse_helper is used in the building of that map
    and free'd afterwards. Not currently used, but may come in handy later on.
  */
  uint32_t* cellse = malloc(4*ncell*sizeof(*cellse));
  uint32_t* cellse_helper = calloc(sizeof(*cellse_helper), ncell);

  x =  malloc(2*nnode*sizeof(*x));
  q =  malloc(4*ncell*sizeof(*q));
  qold = malloc(4*ncell*sizeof(*qold));
  res = malloc(4*ncell*sizeof(*res));
  adt = malloc( ncell*sizeof(*adt));

  for (uint32_t n=0; n<nnode; n++) {
    if (fscanf(fp,"%lf %lf \n",&x[2*n], &x[2*n+1]) != 2) {
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

  /*Populate the cells to edges map*/
  for (uint32_t i = 0; i < nedge; ++i) {
    cellse[4*ecell[2*i] + cellse_helper[ecell[2*i]]] = i;
    cellse_helper[ecell[2*i]]++;
    cellse[4*ecell[2*i+1] + cellse_helper[ecell[2*i+1]]] = i;
    cellse_helper[ecell[2*i+1]]++;
  }
//   for (int i = 0; i < nedge * 2; ++i) {
//     printf("%d:%d\n", i/2, ecell[i]);
//   }
//   printf("-----------------------\n");
/*
   for (int i = 0; i < ncell; ++i) {
      if (cellse_helper[i] != 4) {
        printf("%d - %d\n", i, cellse_helper[i]);
      }
   }
*/
//   free(cellse_helper);
  
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
  for (uint32_t i = 0; i < num_parts; ++i) {
    node_maps[i] = createHashSet(PRIME);
    edge_maps[i] = createHashSet(PRIME);
    cell_maps[i] = createHashSet(PRIME);
    initArr(&ps[i].edges, EDGES_PER_PARTITION);
    initArr(&ps[i].nodes, NODES_PER_PARTITION);
    initArr(&ps[i].cells, CELLS_PER_PARTITION);
    initArr(&ps[i].haloNodes, NODES_PER_PARTITION / 10);
    initArr(&ps[i].haloCells, CELLS_PER_PARTITION / 10);
    initArr(&ps[i].non_halo_cells, CELLS_PER_PARTITION);
  }

  /*Assigning edges to partitions*/
  for (uint32_t i = 0; i < nedge; ++i) {
    uint32_t n1 = npart[edge[2*i]];
//    uint32_t n2 = npart[edge[2*i+1]];
    addToHashSet(edge_maps[n1], i);
//    addToHashSet(edge_maps[n2], i);
  }
  for (uint32_t i = 0; i < num_parts; ++i) {
    ps[i].edges = *toArr(edge_maps[i]);
    destroyHashSet(edge_maps[i]);
  }
  printf("Assigned edges to partitions\n");

  /*
    Assigning cells to partitions, and choosing which ones are halo cells, and non-halo cells.
    Halo cells are cells in which not all nodes belong to the same partition, as determined by the npart
    output of METIS
  */
  
  for (uint32_t i = 0; i < ncell; ++i) {
  //  int n = cpart[i];
//     addToHashSet(cell_maps[n], i);
//    addToArr(&ps[n].cells, i);
    uint32_t p[4];
    for (uint32_t j = 0; j < 4; ++j) {
      p[j] = npart[cell[4*i+j]];
//      addToArr(&ps[p[j]].nodes, cell[4*i+j]);
      addToHashSet(node_maps[p[j]], cell[4*i+j]);
    }
    if (!(p[0] == p[1] && p[0] == p[1] && p[0] == p[2] && p[0] == p[3])) {
      for (uint32_t j = 0; j < 4; ++j) {
        uint32_t t = p[j];
        addToArr(&ps[t].haloCells, i);
        addToHashSet(cell_maps[t], i);
      }
    } else {
      addToArr(&ps[p[0]].non_halo_cells, i);
      addToHashSet(cell_maps[p[0]], i);
    }
  }
 
  /*Clean up the structures populated above and initialise the local cells-to-nodes maps*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    ps[i].nodes = *toArr(node_maps[i]);
    destroyHashSet(node_maps[i]);
    removeDupsArr(&(ps[i].nodes));

    ps[i].cells = *toArr(cell_maps[i]);
    destroyHashSet(cell_maps[i]);
//    removeDupsArr(&(ps[i].cells));
    removeDupsArr(&(ps[i].edges));
    removeDupsArr(&(ps[i].haloCells));
    removeDupsArr(&(ps[i].non_halo_cells));
    ps[i].c2n = malloc(4 * ps[i].cells.len * sizeof(*ps[i].c2n));
  }


  uint32_t p_edges = 0;
  uint32_t p_cells = 0;
  for (uint32_t i = 0; i < num_parts; ++i) {
    /*
      We use hash maps to store the mapping of global cell and node numbers to local numbers.
      We need local cell and node numbers in order to partition the partitions internally.
    */
    ps[i].g2l_nodes = createHashMap(PRIME);
    ps[i].g2l_cells = createHashMap(PRIME);
    ps[i].l2g_nodes = createHashMap(PRIME);
    ps[i].l2g_cells = createHashMap(PRIME);
    uint32_t nodes_added = 0;
    for (uint32_t j = 0; j < ps[i].cells.len; ++j) {
      addToHashMap(ps[i].g2l_cells, ps[i].cells.arr[j], j);
      for (uint32_t k = 0; k < 4; ++k) {
        nodes_added += addToHashMap(ps[i].g2l_nodes, cell[4*ps[i].cells.arr[j] + k], nodes_added);
      }
    }
    for (uint32_t j = 0; j < ps[i].edges.len; ++j) {
      nodes_added += addToHashMap(ps[i].g2l_nodes, edge[2*ps[i].edges.arr[j]], nodes_added);
      nodes_added += addToHashMap(ps[i].g2l_nodes, edge[2*ps[i].edges.arr[j]+1], nodes_added);
    }
    printf("Nodes added %d in partition %d, nodes added in level above: %d\n", nodes_added, i, ps[i].nodes.len);

    /*
      Populate the local cells-to-nodes map, using the hash maps that we filled in above
    */
    for (uint32_t j = 0; j < ps[i].cells.len; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
        uint32_t v = getValue(ps[i].g2l_nodes, cell[4*ps[i].cells.arr[j]+k]);
        ps[i].c2n[4*j + k] = v;
      }
    }

    /*Preparations for partitioning the partition into 2 partitions*/
    uint32_t* pcptr = malloc((ps[i].cells.len + 1) * sizeof(*pcptr));
    for (uint32_t j = 0; j < ps[i].cells.len + 1; ++j) {
      pcptr[j] = 4*j;
    }
    uint32_t* pcpart = malloc(ps[i].cells.len * sizeof(*pcpart));
    uint32_t* pnpart = malloc(nodes_added * sizeof(*pnpart));
    uint32_t nparts = 2;
    METIS_PartMeshNodal((int*)&ps[i].cells.len, (int*)&nodes_added, (int*)pcptr, (int*)ps[i].c2n, NULL, NULL,(int*) &nparts, NULL, NULL, &objval, (int*)pcpart, (int*)pnpart);
    free(pcptr);
    pcptr = NULL;

    /*Here we determine the two internal partitions and the intra-partition halo cells and edges*/
    for (short j = 0; j < 3; ++j) {
      initArr(&ps[i].iparts[j].cells, CELLS_PER_PARTITION / 2);
      initArr(&ps[i].iparts[j].edges, EDGES_PER_PARTITION / 2);
    }
    hash_set* edge_sets[3];
    hash_set* cell_sets[3];
    for (short j = 0; j < 3; ++j) {
      edge_sets[j] = createHashSet(SMALL_PRIME);
      cell_sets[j] = createHashSet(SMALL_PRIME);
    }
    for (uint32_t j = 0; j < ps[i].edges.len; ++j) {
      uint32_t p[2];
    //  printf("ps[%d].edges[%d] = %d\n", i, j, ps[i].edges.arr[j]);
      p[0] = pnpart[getValue(ps[i].g2l_nodes, edge[2*ps[i].edges.arr[j]])];
      p[1] = pnpart[getValue(ps[i].g2l_nodes, edge[2*ps[i].edges.arr[j]+1])];
      if (p[0] != p[1]) {
        addToHashSet(edge_sets[2], ps[i].edges.arr[j]);
      } else {
        addToHashSet(edge_sets[p[0]], ps[i].edges.arr[j]);
      }
    }
    for (short j = 0; j < 3; ++j) {
      ps[i].iparts[j].edges = *toArr(edge_sets[j]);
      destroyHashSet(edge_sets[j]);
    }
    for (uint32_t j = 0; j < ps[i].cells.len; ++j) {
      uint32_t part[4];
      for (uint32_t k = 0; k < 4; ++k) {
        part[k] = pnpart[ps[i].c2n[4*j+k]];
      }
      if (!(part[0] == part[1] && part[0] == part[2] && part[0] == part[3])) {
        addToHashSet(cell_sets[2], j);
//        addToArr(&ps[i].iparts[2].cells, j);
      } else {
        addToHashSet(cell_sets[part[0]], j);
//        addToArr(&ps[i].iparts[part[0]].cells, j);
      }
    }
    for (short j = 0; j < 3; ++j) {
      ps[i].iparts[j].cells = *toArr(cell_sets[j]);
      destroyHashSet(cell_sets[j]);
  //    removeDupsArr(&ps[i].iparts[j].cells);
    }
    p_cells += ps[i].iparts[0].cells.len + ps[i].iparts[1].cells.len + ps[i].iparts[2].cells.len;
    printf("Partition %d is partitioned in partitions of sizes %d, %d and %d intra partition halo cells\n",
           i,
           ps[i].iparts[0].cells.len,
           ps[i].iparts[1].cells.len,
           ps[i].iparts[2].cells.len
          );
    /*
      In this for loop we repeat a similar procedure to above, but we partition each region
      (including the intra-partition halo region into 17 partitions each (TODO: do we need to do this
      for the intra-partition halo as well?)).
      We create an adjacency graph for each of the 3 regions and we colour the bottom-level partitions.
    */
    for (uint32_t j = 0; j < 3; ++j) {
      ipartition* ip = &ps[i].iparts[j];
      ip->g2l_nodes = createHashMap(SMALL_PRIME);
      ip->g2l_cells = createHashMap(SMALL_PRIME);
      ip->l2g_cells = createHashMap(SMALL_PRIME);
      ip->l2g_nodes = createHashMap(SMALL_PRIME);
      ip->c2n = malloc(ip->cells.len * 4 * sizeof(*ip->c2n));
      uint32_t nadded = 0;
      uint32_t nadded_prev = 0;
      for (uint32_t k = 0; k < ip->cells.len; ++k) {
        addToHashMap(ip->g2l_cells, ip->cells.arr[k], k);
        addToHashMap(ip->l2g_cells, k, ip->cells.arr[k]);
        for (uint32_t kk = 0; kk < 4; ++kk) {
          nadded += addToHashMap(ip->g2l_nodes, ps[i].c2n[4*k + kk], nadded);
          if (nadded != nadded_prev) {
            addToHashMap(ip->l2g_nodes, nadded, ps[i].c2n[4*k + kk]);
          }
          nadded_prev = nadded;
        }
      }
      uint32_t na = nadded;

      for (uint32_t k = 0; k < ip->edges.len; ++k) {
        na += addToHashMap(ip->g2l_nodes, getValue(ps[i].g2l_nodes, edge[2*ip->edges.arr[k]]), na);
        na += addToHashMap(ip->g2l_nodes, getValue(ps[i].g2l_nodes, edge[2*ip->edges.arr[k]+1]), na);
      }
      printf("nadded:%d, na - nadded = %d\n", nadded, na - nadded);
      for (uint32_t k = 0; k < ip->cells.len; ++k) {
        for (short kk = 0; kk < 4; ++kk) {
          uint32_t v = getValue(ip->g2l_nodes, ps[i].c2n[4*k + kk]);
          ip->c2n[4*k + kk] = v;
        }
      }
      pcptr = malloc((ip->cells.len + 1) * sizeof(*pcptr));
      for (uint32_t k = 0; k < ip->cells.len + 1; ++k) {
        pcptr[k] = 4*k;
      }
      nparts = BOTTOM_LEVEL_PARTITIONS;
      uint32_t* intcpart = malloc(ip->cells.len * sizeof(*intcpart));
      uint32_t* intnpart = malloc(na * sizeof(*intnpart));
      METIS_PartMeshNodal((int*)&ip->cells.len,(int*)&na, (int*)pcptr, (int*)ip->c2n, NULL, NULL, (int*)&nparts, NULL, NULL, &objval, (int*)intcpart, (int*)intnpart);
      free(pcptr);

      ip->cg = generateGraph(intnpart, ip->c2n, 4, ip->cells.len, nparts);
      colourGraph(ip->cg);
      ip->parts_nodes = malloc(nparts * sizeof(*ip->parts_nodes));
      for (uint32_t k = 0; k < nparts; ++k) {
        initArr(&ip->parts_nodes[k], ip->cells.len / nparts);
      }
      for (uint32_t k = 0; k < nadded; ++k) {
        addToArr(&ip->parts_nodes[intnpart[k]], ps[i].nodes.arr[k]);
      }
      ip->parts_cells = malloc(nparts * sizeof(*ip->parts_cells));
      ip->parts_edges = malloc(nparts * sizeof(*ip->parts_edges));
      hash_set* cell_sets[nparts];
      for (uint32_t k = 0; k < nparts; ++k) {
        cell_sets[k] = createHashSet(SMALL_PRIME);
      }
      for (uint32_t k = 0; k < ip->cells.len; ++k) {
        for (short kk = 0; kk < 4; ++kk) {
          uint32_t n = intnpart[ip->c2n[4 * k + kk]];
          addToHashSet(cell_sets[n], k);
        }
      }
      uint32_t se = 0;
      hash_set* edge_sets[nparts];
      for (uint32_t k = 0; k < nparts; ++k) {
        edge_sets[k] = createHashSet(SMALL_PRIME);
      }
      for (uint32_t k = 0; k < ip->edges.len; ++k) {
        uint32_t p[2];
        p[0] = intnpart[getValue(ip->g2l_nodes, getValue( ps[i].g2l_nodes, edge[ 2 * ip->edges.arr[k] ] ))];
        p[1] = intnpart[getValue(ip->g2l_nodes, getValue( ps[i].g2l_nodes, edge[ 2 * ip->edges.arr[k] + 1 ] ))];
        addToHashSet(edge_sets[p[0]], ip->edges.arr[k]);
      }
      uint32_t bottom_cells = 0;
      for (uint32_t k = 0; k < nparts; ++k) {
        ip->parts_cells[k] = *toArr(cell_sets[k]);
        destroyHashSet(cell_sets[k]);
        
        ip->parts_edges[k] = *toArr(edge_sets[k]);
        destroyHashSet(edge_sets[k]);
        se += ip->parts_edges[k].len;
        printf("bottom level partition %d of partition %d of partition %d has %d cells and %d edges and %d nodes\n", k, j, i, ip->parts_cells[k].len, ip->parts_edges[k].len, ip->parts_nodes[k].len);
        bottom_cells += ip->parts_cells[k].len;
      }
      printf("bottom level cells = %d\n", bottom_cells);
      p_edges += se;
//      printf("internal graph for partition %d of partition %d is:\n", j, i);
//      showGraph(ip->cg);
    }
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
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      initArr(&ps[i].hrCells[j], 64);
    }
  }

  printf("Initialised halo regions\n");

  /*Determine the cells that are shared between every pair of partitions*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    for (uint32_t j = 0; j < ps[i].haloCells.len; ++j) {
      for (uint32_t n = 0; n < ps[i].nneighbours; ++n) {
        if (elemArr(&ps[ps[i].neighbours[n]].haloCells, ps[i].haloCells.arr[j])) {
          addToArr(&ps[i].hrCells[n], ps[i].haloCells.arr[j]);
        }
      }
    }
  }
  /*Clean up the structures from above*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      removeDupsArr(&ps[i].hrCells[j]);
    }
  }

  /*Diagnostic messages, etc...*/
  for (uint32_t i = 0; i < num_parts; ++i) {
    uint32_t total = 0;
    for (uint32_t j = 0; j < ps[i].nneighbours; ++j) {
      printf("Halo region %u-%u has %u cells\n", i, ps[i].neighbours[j], ps[i].hrCells[j].len);
      total += ps[i].hrCells[j].len;
    }
    printf("Total: %u halo cells\n", total);
    printf("----------------------------------\n");
  }

  printf("Colouring partition graph...\n");
  colourGraph(pg);
  printf("Top level graph of mesh partitions:\n");
  showGraph(pg);
  const char* fileName = "meshColoured.dot";
  printf("Writing partition graph to %s ...\n", fileName);
  generateDotGraph(pg, fileName, ps);
  end = clock();
  printf("Edges counted in partitions: %u\n", p_edges);
  printf("Cells counted in partitions: %u\n", p_cells);
  printf("Nodes: %u, Edges: %u, Cells: %u, number of partitions: %u\n", nnode, nedge, ncell, num_parts);
  printf("time taken: %lf seconds\n", (double)(end - start)/CLOCKS_PER_SEC);
  
}
