#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "metis.h"

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

int elem(int*, int, int);
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
    adj_list[i] = malloc(num_parts * sizeof(**adj_list));
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

void toDotColoured(ugraph* g, const char* fileName) {
  if (g->colours == NULL) {
    colourGraph(g);  
  }
  printf("Generating coloured .dot graph in %s\n", fileName);
  FILE* fp = fopen(fileName, "w");
  fprintf(fp, "graph coloured_mesh {\n");
  for (int i = 0; i < g->num_nodes; ++i) {
    fprintf(fp, "\t%d[color=\"%s\"];\n", i, colourNames[g->colours[i]]);
  }
  for (int i = 0; i < g->num_nodes; ++i) {
    for (int j = 0; j < g->adj_sizes[i]; ++j) {
      if (g->adj_list[i][j] > i) {
        fprintf(fp, "\t%d -- %d;\n", i, g->adj_list[i][j]);
      }
    }
  }
  fprintf(fp, "}");
  fclose(fp);
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

int elem(int* arr, int len, int x) {
  for (int i = 0; i < len; ++i) {
    if (arr[i] == x) {
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

void generateDotGraph(int* npart, int* edges, int len, int num_parts, const char* fileName) {
  FILE* fp = fopen(fileName, "w");
  fprintf(fp, "graph mesh {\n");

  int** adj_list =  malloc(num_parts * sizeof(*adj_list));
  int* adj_sizes =  calloc(sizeof(*adj_sizes), num_parts);
  for (int i = 0; i < num_parts; ++i) {
    adj_list[i] = malloc(num_parts * sizeof(**adj_list));
  }
  for (int i = 0; i < len; ++i) {
    int part1 = npart[edges[2*i]];
    int part2 = npart[edges[2*i+1]];
    if (part1 != part2) {
      int pMin = min(part1, part2);
      int pMax = max(part1, part2);
      if (!elem(adj_list[pMin], adj_sizes[pMin], pMax)) {
        adj_list[pMin][adj_sizes[pMin]++] = pMax;
      }
    }
  }
  for (int i = 0; i < num_parts; ++i) {
    for (int j = 0; j < adj_sizes[i]; ++j) {
      fprintf(fp, "\t%d -- %d;\n", i, adj_list[i][j]);
    }
  }
  for (int i = 0; i < num_parts; ++i) {
    free(adj_list[i]);
  }
  free(adj_list);
  free(adj_sizes);
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

typedef struct partition_struct {
  int* edges;
  int* cells;
  int* nodes;
  int* enodes;
  int* ecells;
  int nedges;
  int ncells;
  int nnodes;
  int max_cells;
  int max_edges;
  int max_nodes;
} partition;

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
  printf("Partitioning %d cells in %d partitions\n", ncell, num_parts);
  METIS_PartMeshNodal(&nc, &nn, cptr, cell, NULL, NULL, &num_parts, NULL, NULL, &objval, cpart, npart);

  partition* ps = (partition*)malloc(num_parts * sizeof(*ps));
  for (int i = 0; i < num_parts; ++i) {
    ps[i].nedges = 0;
    ps[i].nnodes = 0;
    ps[i].ncells = 0;
    ps[i].max_edges = EDGES_PER_PARTITION;
    ps[i].max_nodes = NODES_PER_PARTITION;
    ps[i].max_cells = CELLS_PER_PARTITION;
    ps[i].edges = malloc(EDGES_PER_PARTITION * sizeof(*(ps[i].edges)));
    ps[i].nodes = malloc(NODES_PER_PARTITION * sizeof(*(ps[i].nodes)));
    ps[i].cells = malloc(CELLS_PER_PARTITION * sizeof(*(ps[i].cells)));
    ps[i].enodes = malloc(EDGES_PER_PARTITION * 2 * sizeof(*(ps[i].enodes)));
    ps[i].ecells = malloc(EDGES_PER_PARTITION * 2 * sizeof(*(ps[i].ecells)));
  }

  for (int i = 0; i < nedge; ++i) {
    int n1 = cpart[ecell[2*i]];
    int n2 = cpart[ecell[2*i+1]];
    if (ps[n1].nedges >= ps[n1].max_edges - 1) {
      ps[n1].max_edges += 128;
      ps[n1].edges = realloc(ps[n1].edges, ps[n1].max_edges * sizeof(*(ps[n1].edges)));
    }
    ps[n1].edges[ps[n1].nedges++] = i;

    if (ps[n2].nedges >= ps[n2].max_edges - 1) {
      ps[n2].max_edges += 128;
      ps[n2].edges = realloc(ps[n2].edges, ps[n2].max_edges * sizeof(*(ps[n2].edges)));
    }
    ps[n2].edges[ps[n2].nedges++] = i;
  }
  printf("Assigned edges to partitions\n");

  for (int i = 0; i < ncell; ++i) {
    int n = cpart[i]; // n is partition number
    if (ps[n].ncells >= ps[n].max_cells - 1) {
      ps[n].max_cells += 128;
      ps[n].cells = realloc(ps[n].cells, ps[n].max_cells * sizeof(*ps[n].cells));
    }
    ps[n].cells[ps[n].ncells++] = i;
    
    for (int n = 0; n < 4; ++n) {
      int np = npart[cell[4*i+n]];
      if (ps[np].nnodes == ps[np].max_nodes - 1) {
        ps[np].max_nodes += 128;
        ps[np].nodes = realloc(ps[np].nodes, ps[np].max_nodes * sizeof(*ps[np].nodes));
      }
      ps[np].nodes[ps[np].nnodes++] = cell[4*i];
    }
  }
  
  for (int i = 0; i < num_parts; ++i) {
    ps[i].nnodes = removeDups(ps[i].nodes, ps[i].nnodes);
    ps[i].nodes = realloc(ps[i].nodes, ps[i].nnodes * sizeof(*ps[i].nodes));
    
    ps[i].ncells = removeDups(ps[i].cells, ps[i].ncells);
    ps[i].cells = realloc(ps[i].cells, ps[i].ncells * sizeof(*ps[i].cells));
    
    ps[i].nedges = removeDups(ps[i].edges, ps[i].nedges);
    ps[i].edges = realloc(ps[i].edges, ps[i].nedges * sizeof(*ps[i].edges));
  
  }

  int part_nodes = 0;
  int part_cells = 0;
  int part_edges = 0;
  for (int i = 0; i < num_parts; ++i) {
    part_nodes += ps[i].nnodes;
    part_cells += ps[i].ncells;
    part_edges += ps[i].nedges;
    printf("partition %d has: %d edges, %d nodes, %d cells\n", i, ps[i].nedges, ps[i].nnodes, ps[i].ncells);
  }
//  ugraph* partitionGraph = generateGraph(npart, edge, nedge*2, num_parts);
  printf("Generating partition graph\n");
  ugraph* partitionGraph = generateGraph(npart, cell, 4, ncell, num_parts);
  if (!partitionGraph) {
    printf("ERROR! partition graph is NULL!\n");
    return 1;
  }
  showGraph(partitionGraph);
  printf("Colouring partition graph\n");
  colourGraph(partitionGraph);
  toDotColoured(partitionGraph, "meshColoured.dot");
//  generateDotGraph(npart, edge, nedge*2, num_parts, "mesh.dot");
  end = clock();
  printf("Nodes: %d, Edges: %d, Cells: %d\n", nnode, nedge, ncell);
  printf("Halo data:\n Nodes: %d, Edges: %d, Cells: %d\n", part_nodes - nnode, part_edges - nedge, part_cells - ncell);
  printf("time taken: %lf seconds\n", (double)(end - start)/CLOCKS_PER_SEC);
  
}
