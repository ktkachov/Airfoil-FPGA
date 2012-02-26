#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "metis.h"
#include "hash_map.h"

#define EDGES_PER_PARTITION 80000
#define NODES_PER_PARTITION (EDGES_PER_PARTITION / 2)
#define CELLS_PER_PARTITION (EDGES_PER_PARTITION / 2)

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

typedef struct graph_struct {
  int** adj_list;
  int* adj_sizes;
  int* colours;
  int num_nodes;
} ugraph;

typedef struct part_array_struct {
  int* arr;
  int len;
  int maxLen;
} arr_t;

typedef struct partition_struct {
  int* enodes;
  int* ecells;

  int* neighbours;
  int nneighbours;

  int* c2n; /*Cells to nodes map*/

  hash_map* node_ind; /*mapping of local node number to global node number*/

  arr_t haloCells;
  arr_t haloNodes;
  arr_t edges;
  arr_t cells;
  arr_t nodes;
  arr_t* hrCells; /*Halo region cells*/
} partition;

int elem(int*, int, int);
int elemArr(arr_t*, int);
void colourGraph(ugraph*);
int removeDups(int*, int);
int min(int, int);
int max(int, int);
void showGraph(ugraph*);

ugraph* generateGraph(int* npart, int* map, int dim, int len, int num_parts) {
  int** adj_list = malloc(num_parts * sizeof(*adj_list));
  if (!adj_list) {
    fprintf(stderr, "ERROR! could not allocate memory for adjacency list\n");
  }
  int* adj_sizes = calloc(sizeof(*adj_sizes), num_parts);
  for (int i = 0; i < num_parts; ++i) {
    adj_list[i] = malloc(num_parts * sizeof(*adj_list[i]));
  } 
  for (int i = 0; i < len; ++i) {
    int p[dim];
    for (int j = 0; j < dim; ++j) {
      p[j] = npart[map[dim*i+j]];
    }
    for (int j = 0; j < dim; ++j) {
      for (int k = 0; k < dim; ++k) {
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
  int nv = graph->num_nodes; //Number of vertices, also maximum colours
  int* colours = malloc(graph->num_nodes * sizeof(*colours));
  for (int i = 0; i < nv; ++i) {
    colours[i] = -1;
  }
  int ncolours = 1;
  for (int i = 0; i < nv; ++i) {
    int c = 0;
    int repeat = 1;
    while (c < ncolours && repeat) {
      int a = 1; //colour availble
      for (int v = 0; v < graph->adj_sizes[i] && a; ++v) {
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

inline void print_array(double* data, int len, const char* file_name) {
  FILE* flog;
  flog = fopen(file_name, "w");
  for (int i = 0; i < len; ++i) {
    fprintf(flog, "%lf\n", data[i]);
  }
  fclose(flog);
  
}

int removeDups(int* arr, int len) {
  int i, j;
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
  int i, j;
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

inline int elem(int* arr, int len, int x) {
  for (int i = 0; i < len; ++i) {
    if (arr[i] == x) {
      return 1;
    }
  }
  return 0;
}

inline int elemArr(arr_t* a, int x) {
  for (int i = 0; i < a->len; ++i) {
    if (a->arr[i] == x) {
      return 1;
    }
  }
  return 0;
}

inline int min(int a, int b) {
  return a < b ? a : b;
}

inline int max(int a, int b) {
  return a > b ? a : b;
}

void generateDotGraph(ugraph* g, const char* fileName, partition* ps) {
  FILE* fp = fopen(fileName, "w");
  fprintf(fp, "graph mesh {\n");

  if (g->colours != NULL) {
    for (int i = 0; i < g->num_nodes; ++i) {
      fprintf(fp, "\t%d[color=\"%s\"];\n", i, colourNames[g->colours[i]]);
    }
  }  
  for (int i = 0; i < g->num_nodes; ++i) {
    for (int j = 0; j < g->adj_sizes[i]; ++j) {
      if (g->adj_list[i][j] > i) {
        fprintf(fp, "\t%d -- %d [label=%d];\n", i, g->adj_list[i][j], ps[i].hrCells[j].len);
      }
    }
  }
  fprintf(fp, "}");
  fclose(fp);
}


void showGraph(ugraph* g) {
  for (int i = 0; i < g->num_nodes; ++i) {
    printf("%d : ", i);
    for (int j = 0 ; j < g->adj_sizes[i]; ++j) {
      printf("%d, ", g->adj_list[i][j]);
    }
    printf("\n");
  }
}


void addToArr(arr_t* a, int e) {
  if (a->len == a->maxLen - 1) {
    a->maxLen += 64;
    a->arr = realloc(a->arr, a->maxLen * sizeof(*(a->arr)));
  }
  a->arr[a->len++] = e;
}

void initArr(arr_t* a, int size) {
  a->len = 0;
  a->maxLen = size;
  a->arr = malloc(a->maxLen * sizeof(*(a->arr)));
}


int main(int argc, char* argv[]) {

  int *becell, *ecell, *bound, *bedge, *edge, *cell;
  double *x, *q, *qold, *adt, *res;

  int nnode,ncell,nedge,nbedge;
  
  // read in grid

  printf("reading in grid \n");
  FILE *fp;
  if ( (fp = fopen("./new_grid.dat","r")) == NULL) { ///new_grid.dat
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell = malloc(4*ncell*sizeof(*cell));
  edge = malloc(2*nedge*sizeof(*edge));
  ecell = malloc(2*nedge*sizeof(*ecell));
  bedge =  malloc(2*nbedge*sizeof(*bedge));
  becell = malloc( nbedge*sizeof(*becell));
  bound =  malloc( nbedge*sizeof(*bound));

  int* cellse = malloc(4*ncell*sizeof(*cellse));
  int* cellse_helper = calloc(sizeof(*cellse_helper), ncell);

  x =  malloc(2*nnode*sizeof(*x));
  q =  malloc(4*ncell*sizeof(*q));
  qold = malloc(4*ncell*sizeof(*qold));
  res = malloc(4*ncell*sizeof(*res));
  adt = malloc( ncell*sizeof(*adt));

  for (int n=0; n<nnode; n++) {
    if (fscanf(fp,"%lf %lf \n",&x[2*n], &x[2*n+1]) != 2) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<ncell; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&cell[4*n ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<nedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (int n=0; n<nbedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&bedge[2*n],&bedge[2*n+1],
                                   &becell[n], &bound[n]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  fclose(fp);


/*Populate the cells to edges map*/
  for (int i = 0; i < nedge; ++i) {
    cellse[4*ecell[2*i] + cellse_helper[ecell[2*i]]++] = i;
    cellse[4*ecell[2*i+1] + cellse_helper[ecell[2*i+1]]++] = i;
  }
  free(cellse_helper);
  /*
  for (int i = 0; i < ncell * 4 ; i+=4) {
    printf("%d, %d, %d, %d\n", cellse[i], cellse[i+1], cellse[i+2], cellse[i+3]);
  }
*/
  clock_t start, end;
  start = clock();
 
  //Cell pointer for METIS
  int* cptr = malloc((ncell+1) * sizeof(*cptr));
  for (int i = 0; i < ncell+1; ++i) {
    cptr[i] =4*i;
  }
  int objval;
  int *cpart = malloc(ncell * sizeof(*cpart));
  int *npart = malloc(nnode * sizeof(*npart));
  int num_parts = ncell / CELLS_PER_PARTITION;
  int nn = nnode;
  int nc = ncell;
  printf("Partitioning %d cells in %d partitions...\n", ncell, num_parts);
  METIS_PartMeshNodal(&nc, &nn, cptr, cell, NULL, NULL, &num_parts, NULL, NULL, &objval, cpart, npart);

  partition* ps = (partition*)malloc(num_parts * sizeof(*ps));
  for (int i = 0; i < num_parts; ++i) {

    initArr(&ps[i].edges, EDGES_PER_PARTITION);
    initArr(&ps[i].nodes, NODES_PER_PARTITION);
    initArr(&ps[i].cells, CELLS_PER_PARTITION);
    initArr(&ps[i].haloNodes, NODES_PER_PARTITION / 10);
    initArr(&ps[i].haloCells, CELLS_PER_PARTITION / 10);

  }

  for (int i = 0; i < nedge; ++i) {
    int n1 = cpart[ecell[2*i]];
    int n2 = cpart[ecell[2*i+1]];
    addToArr(&ps[n1].edges, i);
    addToArr(&ps[n2].edges, i);
  }
  printf("Assigned edges to partitions\n");

  for (int i = 0; i < ncell; ++i) {
    int n = cpart[i]; // n is partition number
    addToArr(&ps[n].cells, i);
    
    int p[4];
    for (int j = 0; j < 4; ++j) {
      p[j] = npart[cell[4*i+j]];
    }
    if (!(p[0] == p[1] && p[0] == p[1] && p[0] == p[2] && p[0] == p[3])) {
      for (int j = 0; j < 4; ++j) {
        int t = p[j];
        addToArr(&ps[t].haloCells, i);
      }
    }
    for (int n = 0; n < 4; ++n) {
      int np = npart[cell[4*i+n]];
      addToArr(&ps[np].nodes, cell[4*i]);
    }
  }
 
  for (int i = 0; i < num_parts; ++i) {
    removeDupsArr(&ps[i].nodes);
    removeDupsArr(&ps[i].cells);
    removeDupsArr(&ps[i].edges);
    removeDupsArr(&ps[i].haloCells);
    ps[i].c2n = malloc(4 * ps[i].cells.len * sizeof(*ps[i].c2n));
  }

  printf("Assigning cells to nodes maps to partitions...\n");
  for (int i = 0; i < num_parts; ++i) {
    ps[i].node_ind = createHashMap(ps[i].nodes.len);
    for (int j = 0; j < ps[i].cells.len; ++j) {
      for (int k = 0; k < 4; ++k) {
        ps[i].c2n[4*j+k] = cell[4*ps[i].cells.arr[j] + k];
      }
    }
  }

  for (int i = 0; i < num_parts; ++i) {
    ps[i].node_ind = createHashMap(ps[i].nodes.len);
    for (int j = 0; j < ps[i].nodes.len; ++j) {
      addToHashMap(ps[i].node_ind, ps[i].nodes.arr[j], j);
    }
  }

  int part_nodes = 0;
  int part_cells = 0;
  int part_edges = 0;
  for (int i = 0; i < num_parts; ++i) {
    part_nodes += ps[i].nodes.len;
    part_cells += ps[i].cells.len;
    part_edges += ps[i].edges.len;
    printf("partition %d has: %d edges, %d nodes, %d cells, %d halo cells\n", i, ps[i].edges.len, ps[i].nodes.len, ps[i].cells.len, ps[i].haloCells.len);
  }
//  ugraph* partitionGraph = generateGraph(npart, edge, nedge*2, num_parts);
  printf("Generating partition graph\n");
  ugraph* pg = generateGraph(npart, cell, 4, ncell, num_parts);
  if (!pg) {
    printf("ERROR! partition graph is NULL!\n");
    return 1;
  }

  for (int i = 0; i < num_parts; ++i) {
    ps[i].neighbours = malloc(pg->adj_sizes[i] * sizeof(*(ps[i].neighbours)));
    ps[i].nneighbours = pg->adj_sizes[i];
    memcpy(ps[i].neighbours, pg->adj_list[i], pg->adj_sizes[i] * sizeof(*(ps[i].neighbours)));
    ps[i].hrCells = malloc(ps[i].nneighbours * sizeof(*ps[i].hrCells));
    for (int j = 0; j < ps[i].nneighbours; ++j) {
      initArr(&ps[i].hrCells[j], 64);
    }
  }

  printf("Initialised halo regions\n");

  for (int i = 0; i < num_parts; ++i) {
    for (int j = 0; j < ps[i].haloCells.len; ++j) {
      for (int n = 0; n < ps[i].nneighbours; ++n) {
        if (elemArr(&ps[ps[i].neighbours[n]].haloCells, ps[i].haloCells.arr[j])) {
          addToArr(&ps[i].hrCells[n], ps[i].haloCells.arr[j]);
        }
      }
    }
  }

  for (int i = 0; i < num_parts; ++i) {
    for (int j = 0; j < ps[i].nneighbours; ++j) {
      removeDupsArr(&ps[i].hrCells[j]);
    }
  }

  for (int i = 0; i < num_parts; ++i) {
    int total = 0;
    for (int j = 0; j < ps[i].nneighbours; ++j) {
      printf("Halo region %d-%d has %d cells\n", i, ps[i].neighbours[j], ps[i].hrCells[j].len);
      total += ps[i].hrCells[j].len;
    }
    printf("Total: %d halo cells\n", total);
    printf("----------------------------------\n");
  }

  showGraph(pg);
  printf("Colouring partition graph...\n");
  colourGraph(pg);
  const char* fileName = "meshColoured.dot";
  printf("Writing partition graph to %s ...\n", fileName);
  generateDotGraph(pg, fileName, ps);
//  generateDotGraph(npart, edge, nedge*2, num_parts, "mesh.dot");
  end = clock();
  printf("Nodes: %d, Edges: %d, Cells: %d\n", nnode, nedge, ncell);
  printf("time taken: %lf seconds\n", (double)(end - start)/CLOCKS_PER_SEC);
  
}
