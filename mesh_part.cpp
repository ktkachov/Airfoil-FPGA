#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "metis.h"

#define EDGES_PER_PARTITION 80000
#define NODES_PER_PARTITION 40000
#define CELLS_PER_PARTITION 40000

double gam;
double gm1;
double cfl;
double eps;
double qinf[4];

const char* colourNames[20] = {
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

ugraph* generateGraph(int* npart, int* edges, int len, int num_parts) {
  
  int** adj_list = (int**) malloc(num_parts * sizeof(int));
  int* adj_sizes = (int*) calloc(sizeof(int), num_parts);
  for (int i = 0; i < num_parts; ++i) {
    adj_list[i] = (int*) malloc(num_parts * sizeof(num_parts));
  }
  for (int i = 0; i < len; ++i) {
    int part1 = npart[edges[2*i]];
    int part2 = npart[edges[2*i+1]];
    if (part1 != part2) {
      if (!elem(adj_list[part1], adj_sizes[part1], part2)) {
        adj_list[part1][adj_sizes[part1]++] = part2;
      }
      if (!elem(adj_list[part2], adj_sizes[part2], part1)) {
        adj_list[part2][adj_sizes[part2]++] = part1;
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
    printf("Node %d has colour %d\n", i, g->colours[i]);
    fprintf(fp, "%d[color=\"%s\"];\n", i, colourNames[g->colours[i]]);
  }
  for (int i = 0; i < g->num_nodes; ++i) {
    for (int j = 0; j < g->adj_sizes[i]; ++j) {
      if (g->adj_list[i][j] > i) {
        fprintf(fp, "%d -- %d;\n", i, g->adj_list[i][j]);
      }
    }
  }
  fprintf(fp, "}");
  fclose(fp);
}

void colourGraph(ugraph* graph) {
  int nv = graph->num_nodes; //Number of vertices, also maximum colours
  int* colours = (int*)malloc(graph->num_nodes * sizeof(int));
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
  for (i = 1; i < len; i++) {
    if (arr[i] != arr[j]) {
      j++;
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

  int** adj_list = (int**) malloc(num_parts * sizeof(int));
  int* adj_sizes = (int*) calloc(sizeof(int), num_parts);
  for (int i = 0; i < num_parts; ++i) {
    adj_list[i] = (int*) malloc(num_parts * sizeof(num_parts));
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
      fprintf(fp, "%d -- %d;\n", i, adj_list[i][j]);
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

typedef struct partition_struct {
  int* edges;
  int* cells;
  int* nodes;
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

  int nnode,ncell,nedge,nbedge,niter;
  double rms;
  
  //timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;

  // read in grid

  printf("reading in grid \n");

  FILE *fp;
  if ( (fp = fopen("./new_grid.dat","r")) == NULL) { ///new_grid.dat
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell = (int *) malloc(4*ncell*sizeof(int));
  edge = (int *) malloc(2*nedge*sizeof(int));
  ecell = (int *) malloc(2*nedge*sizeof(int));
  bedge = (int *) malloc(2*nbedge*sizeof(int));
  becell = (int *) malloc( nbedge*sizeof(int));
  bound = (int *) malloc( nbedge*sizeof(int));

  int* cellse = (int*) malloc(4*ncell*sizeof(int));
  int* cellse_helper = (int*)calloc(sizeof(int), ncell);

  x = (double *) malloc(2*nnode*sizeof(double));
  q = (double *) malloc(4*ncell*sizeof(double));
  qold = (double *) malloc(4*ncell*sizeof(double));
  res = (double *) malloc(4*ncell*sizeof(double));
  adt = (double *) malloc( ncell*sizeof(double));

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

  // set constants and initialise flow field and residual
  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;


  double mach = 0.4f;
  double alpha = 3.0f*atan(1.0f)/45.0f;
  double p = 1.0f;
  double r = 1.0f;
  double u = sqrt(gam*p/r)*mach;
  double e = p/(r*gm1) + 0.5f*u*u;

  qinf[0] = r;
  qinf[1] = r*u;
  qinf[2] = 0.0f;
  qinf[3] = r*e;

  for (int n=0; n<ncell; n++) {
    for (int m=0; m<4; m++) {
        q[4*n+m] = qinf[m];
      res[4*n+m] = 0.0f;
      
    }
  }

  niter = 1000;
//   printf("Before giant for loop\n");

/*Populate the cells to edges map*/
  for (int i = 0; i < nedge; ++i) {
    cellse[4*ecell[2*i] + cellse_helper[ecell[2*i]]++] = i;
    cellse[4*ecell[2*i+1] + cellse_helper[ecell[2*i+1]]++] = i+1;
  }
  free(cellse_helper);
  /*
  for (int i = 0; i < ncell * 4 ; i+=4) {
    printf("%d, %d, %d, %d\n", cellse[i], cellse[i+1], cellse[i+2], cellse[i+3]);
  }
*/
  clock_t start, end;
  start = clock();
 
  int* eptr = (int*)malloc((nedge+1) * sizeof(int));
  for (int i = 0; i < nedge+1; ++i) {
    eptr[i] =2*i;
  }
  int objval;
  int *epart = (int*)malloc(nedge * sizeof(int));
  int *npart = (int*)malloc(nnode * sizeof(int));
  int num_parts = nedge / EDGES_PER_PARTITION;
  //idx_t options[METIS_NOPTIONS];
  //options[METIS_OPTION_NUMBERING] = 0;
  int nn = nnode;
  int ne = nedge;
  printf("About to start partitioning %d edges in %d partitions\n", nedge, num_parts);
  METIS_PartMeshNodal(&ne, &nn, eptr, edge, NULL, NULL, &num_parts, NULL, /*options*/NULL, &objval, epart, npart);

  partition* ps = (partition*)malloc(num_parts * sizeof(partition));
  for (int i = 0; i < num_parts; ++i) {
    ps[i].nedges = 0;
    ps[i].nnodes  = 0;
    ps[i].ncells = 0;
    ps[i].max_edges = EDGES_PER_PARTITION;
    ps[i].max_nodes = NODES_PER_PARTITION;
    ps[i].max_cells = CELLS_PER_PARTITION;
    ps[i].edges = (int*)malloc(EDGES_PER_PARTITION * sizeof(int));
    ps[i].nodes = (int*)malloc(NODES_PER_PARTITION * sizeof(int));
    ps[i].cells = (int*)malloc(CELLS_PER_PARTITION * sizeof(int));
  }
  for (int i = 0; i < ne; ++i) {
    int n = epart[i]; // n is partition number
    if(ps[n].nedges == ps[n].max_edges) {
      ps[n].max_edges += 64;
      ps[n].edges = (int*)realloc(ps[n].edges, ps[n].max_edges * sizeof(int));
    }
    ps[n].edges[ps[n].nedges++] = i;
  }
  
  printf("Assigned edges to partitions\n");

  for (int i = 0; i < nn; ++i) {
    int n = npart[i];
    if(ps[n].nnodes == ps[n].max_nodes) {
      ps[n].max_nodes += 64;
      ps[n].nodes = (int*)realloc(ps[n].nodes, ps[n].max_nodes * sizeof(int));
    }
    ps[n].nodes[ps[n].nnodes++] = i;
  }
    

  printf("Assigned nodes to partitions\n");

  for (int i = 0; i < ncell; ++i) {
    for (int j = 0; j < 4; ++j) {
      int edge = cellse[4*i+j];
      if (edge > nedge) {
        printf("ERROR! Edge %d\n", edge);
      }
      int part = epart[edge];
      if (ps[part].ncells == ps[part].max_cells) {
        ps[part].max_cells += 64;
        ps[part].cells = (int*)realloc(ps[part].cells, ps[part].max_cells * sizeof(int));
      }
      ps[part].cells[ps[part].ncells++] = i;
    }
  }
/*
  for (int i = 0; i < nedge*2; i+=2) {
    int n = epart[i/2];
    if (ps[n].ncells > ps[n].max_cells-2) {
      ps[n].max_cells += 64;
      ps[n].cells = (int*)realloc(ps[n].cells, ps[n].max_cells * sizeof(int));
    }
    ps[n].cells[ps[n].ncells++] = ecell[i];
    ps[n].cells[ps[n].ncells++] = ecell[i+1];
  }
*/
  //Clean up and remove duplicates
  for (int i = 0; i < num_parts; ++i) {
    printf("Nodes before in partition %d: %d\n", i, ps[i].nnodes);
    ps[i].nnodes = removeDups(ps[i].nodes, ps[i].nnodes);
    ps[i].nodes = (int*)realloc(ps[i].nodes, ps[i].nnodes * sizeof(int));
    printf("Nodes after in partition %d: %d\n", i, ps[i].nnodes);

    printf("Cells before in partition %d: %d\n",i,ps[i].ncells);
    ps[i].ncells = removeDups(ps[i].cells, ps[i].ncells);
    ps[i].cells = (int*)realloc(ps[i].cells, ps[i].ncells * sizeof(int));
    printf("Cells after in partition %d: %d\n", i, ps[i].ncells);
  }

  printf("Assigned cells to partitions\n");


  int part_nodes = 0;
  int part_cells = 0;
  int part_edges = 0;
  for (int i = 0; i < num_parts; ++i) {
    part_nodes += ps[i].nnodes;
    part_cells += ps[i].ncells;
    part_edges += ps[i].nedges;
    printf("partition %d has: %d edges, %d nodes, %d cells\n", i, ps[i].nedges, ps[i].nnodes, ps[i].ncells);
  }
/*
  for (int i = 0; i < num_parts; ++i) {
    printf("printing partition %d of %d\n", i, num_parts);
    for (int j = 0; j < ps[i].nedges; ++j) {
      printf("%d\n", ps[i].edges[j]);
    }
    printf("----------------------------------\n");
  }
*/
/*
  for (int i = 0; i < nedge; ++i) {
    printf("%d\n", epart[i]);
  }
 */
/*
  for (int i = 0; i < nnode; ++i) {
    printf("%d\n", npart[i]);
  }
*/
  ugraph* partitionGraph = generateGraph(npart, edge, nedge*2, num_parts);
  if (partitionGraph == NULL) {
    printf("ERROR! partition graph is NULL!\n");
    return 1;
  }
  colourGraph(partitionGraph);
  toDotColoured(partitionGraph, "meshColoured.dot");
//  generateDotGraph(npart, edge, nedge*2, num_parts, "mesh.dot");
  end = clock();
  printf("Nodes: %d, Edges: %d, Cells: %d\n", nnode, nedge, ncell);
  printf("Halo data:\n Nodes: %d, Edges: %d, Cells: %d\n", part_nodes - nnode, part_edges - nedge, part_cells - ncell);
  printf("time taken: %lf seconds\n", (double)(end - start)/CLOCKS_PER_SEC);
  
}
