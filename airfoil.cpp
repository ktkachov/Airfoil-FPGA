#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "airfoil_kernels.h"

float gam;
float gm1;
float cfl;
float eps;
float qinf[4];

inline void print_array(float* data, int len, const char* file_name) {
  FILE* flog;
  flog = fopen(file_name, "w");
  for (int i = 0; i < len; ++i) {
    fprintf(flog, "%lf\n", data[i]);
  }
  fclose(flog);
  
}

int main(int argc, char **argv){

  int *becell, *ecell, *bound, *bedge, *edge, *cell;
  float *x, *q, *qold, *adt, *res;

  int nnode,ncell,nedge,nbedge,niter;
  float rms;
  
  // read in grid

  printf("reading in grid \n");

  FILE *fp;
  if ( (fp = fopen("./new_grid.dat","r")) == NULL) { ///new_grid.dat
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  cell = new int[4*ncell];
  edge = new int[2*nedge];
  ecell = new int[2*nedge];
  bedge = new int[2*nbedge];
  becell = new int[nbedge];
  bound = new int[nbedge];

  x = new float[2*nnode];
  q = new float[4*ncell];  
  qold = new float[4*ncell];
  res = new float[4*ncell];
  adt = new float[ncell];

  for (int n=0; n<nnode; n++) {
    if (fscanf(fp,"%f %f \n",&x[2*n], &x[2*n+1]) != 2) {
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

  printf("initialising flow field \n");

  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;


  float mach = 0.4f;
  //float alpha = 3.0f*atan(1.0f)/45.0f;
  float p = 1.0f;
  float r = 1.0f;
  float u = sqrt(gam*p/r)*mach;
  float e = p/(r*gm1) + 0.5f*u*u;

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

  clock_t start, end;
  start = clock();
 
  double times[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
 
  for(int iter=1; iter<=niter; iter++) {

// save old flow solution
    clock_t before_clock = clock();
    for (int i = 0; i < ncell; ++i) {
      save_soln(&q[4*i], &qold[4*i]);
    }
    times[0] += clock() - before_clock;
// predictor/corrector update loop


    for(int k=0; k<2; k++) {

// calculate area/timstep
      before_clock = clock();
      for (int i = 0; i < ncell; ++i) {
        adt_calc(&x[2*cell[4*i]], &x[2*cell[4*i+1]], &x[2*cell[4*i+2]], &x[2*cell[4*i+3]], &q[4*i], &adt[i]);
      }
      times[1] += clock() - before_clock;

// calculate flux residual

      before_clock = clock();
      for (int i = 0; i < nedge; ++i) {
        res_calc(&x[2*edge[2*i]], &x[2*edge[2*i+1]], &q[4*ecell[2*i]], &q[4*ecell[2*i+1]], &adt[ecell[2*i]], &adt[ecell[2*i+1]], &res[4*ecell[2*i]], &res[4*ecell[2*i+1]]);
      }
      times[2] += clock() - before_clock;

      before_clock = clock();
      for (int i = 0; i < nbedge; ++i) {
        bres_calc(&x[2*bedge[2*i]], &x[2*bedge[2*i+1]], &q[4*becell[i]], &adt[becell[i]], &res[4*becell[i]], &bound[i]);
      }
      times[3] += clock() - before_clock;


// update flow field

      rms = 0.0;

      before_clock = clock();
      for (int i = 0; i < ncell; ++i) {
        update(&qold[4*i], &q[4*i], &res[4*i], &adt[i], &rms);
      }
      times[4] += clock() - before_clock;
    }

// print iteration history

    rms = sqrt(rms/(float) ncell);

    if (iter%100 == 0)
      printf(" %d %10.5e \n",iter,rms);
  }
  
  end = clock();
  printf("time taken: %lf seconds\n", (float)(end - start)/CLOCKS_PER_SEC);
  printf("save_soln: %lf seconds\n", (float)times[0]/CLOCKS_PER_SEC);
  printf("adt_calc: %lf seconds\n", (float)times[1]/CLOCKS_PER_SEC);
  printf("res_calc: %lf seconds\n", (float)times[2]/CLOCKS_PER_SEC);
  printf("bres_calc: %lf seconds\n", (float)times[3]/CLOCKS_PER_SEC);
  printf("update: %lf seconds\n", (float)times[4]/CLOCKS_PER_SEC);
}
