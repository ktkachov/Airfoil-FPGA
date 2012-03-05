#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

double gam;
double gm1;
double cfl;
double eps;
double qinf[4];

inline void print_array(double* data, long len, const char* file_name) {
  FILE* flog;
  flog = fopen(file_name, "w");
  for (long i = 0; i < len; ++i) {
    fprintf(flog, "%lf\n", data[i]);
  }
  fclose(flog);
  
}

inline void adt_calc(double *x1,double *x2,double *x3,double *x4,double *q,double *adt){
  double dx,dy, ri,u,v,c;

  ri = 1.0f/q[0];
//   printf("ri=%lf, q[0]=%lf\n", ri, q[0]);
  u = ri*q[1];
  v = ri*q[2];
  c = sqrt(gam*gm1*(ri*q[3]-0.5f*(u*u+v*v)));

  dx = x2[0] - x1[0];
  dy = x2[1] - x1[1];
  *adt = fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x3[0] - x2[0];
  dy = x3[1] - x2[1];
  *adt += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x4[0] - x3[0];
  dy = x4[1] - x3[1];
  *adt += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x1[0] - x4[0];
  dy = x1[1] - x4[1];
  *adt += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  *adt = (*adt) / cfl;
}

inline void bres_calc(double *x1, double *x2, double *q1,
                      double *adt1,double *res1,long *bound) {
  double dx,dy,mu, ri, p1,vol1, p2,vol2, f;

  dx = x1[0] - x2[0];
  dy = x1[1] - x2[1];

  ri = 1.0f/q1[0];
  p1 = gm1*(q1[3]-0.5f*ri*(q1[1]*q1[1]+q1[2]*q1[2]));

  if (*bound==1) {
    res1[1] += + p1*dy;
    res1[2] += - p1*dx;
  }
  else {
    vol1 = ri*(q1[1]*dy - q1[2]*dx);

    ri = 1.0f/qinf[0];
    p2 = gm1*(qinf[3]-0.5f*ri*(qinf[1]*qinf[1]+qinf[2]*qinf[2]));
    vol2 = ri*(qinf[1]*dy - qinf[2]*dx);

    mu = (*adt1)*eps;

    f = 0.5f*(vol1* q1[0] + vol2* qinf[0] ) + mu*(q1[0]-qinf[0]);
    res1[0] += f;
    f = 0.5f*(vol1* q1[1] + p1*dy + vol2* qinf[1] + p2*dy) + mu*(q1[1]-qinf[1]);
    res1[1] += f;
    f = 0.5f*(vol1* q1[2] - p1*dx + vol2* qinf[2] - p2*dx) + mu*(q1[2]-qinf[2]);
    res1[2] += f;
    f = 0.5f*(vol1*(q1[3]+p1) + vol2*(qinf[3]+p2) ) + mu*(q1[3]-qinf[3]);
    res1[3] += f;
  }
}

inline void res_calc(double *x1, double *x2, double *q1, double *q2,
                     double *adt1,double *adt2,double *res1,double *res2) {
  
//    printf("adt1 = %lf, adt2= %lf, res1[0] = %lf, res2[0] = %lf\n", *adt1, *adt2, res1[0], res2[0]);
//   
//    for (int i = 0; i < 4; ++i) {
//       printf("q1[%d] = %lf, q2[%d] = %lf\n", i, q1[i], i, q2[i]);
//    }
   
  double dx,dy,mu, ri, p1,vol1, p2,vol2, f;

  dx = x1[0] - x2[0];
  dy = x1[1] - x2[1];

  ri = 1.0f/q1[0];
  p1 = gm1*(q1[3]-0.5f*ri*(q1[1]*q1[1]+q1[2]*q1[2]));
  vol1 = ri*(q1[1]*dy - q1[2]*dx);

  ri = 1.0f/q2[0];
  p2 = gm1*(q2[3]-0.5f*ri*(q2[1]*q2[1]+q2[2]*q2[2]));
  vol2 = ri*(q2[1]*dy - q2[2]*dx);

  mu = 0.5f*((*adt1)+(*adt2))*eps;

  f = 0.5f*(vol1* q1[0] + vol2* q2[0] ) + mu*(q1[0]-q2[0]);
  res1[0] += f;
  res2[0] -= f;
  f = 0.5f*(vol1* q1[1] + p1*dy + vol2* q2[1] + p2*dy) + mu*(q1[1]-q2[1]);
  res1[1] += f;
  res2[1] -= f;
  f = 0.5f*(vol1* q1[2] - p1*dx + vol2* q2[2] - p2*dx) + mu*(q1[2]-q2[2]);
  res1[2] += f;
  res2[2] -= f;
  f = 0.5f*(vol1*(q1[3]+p1) + vol2*(q2[3]+p2) ) + mu*(q1[3]-q2[3]);
  res1[3] += f;
  res2[3] -= f;
}

inline void save_soln(double *q, double *qold){
  for (int n=0; n<4; n++) qold[n] = q[n];
}

inline void update(double *qold, double *q, double *res, double *adt, double *rms){
  double del, adti;

  adti = 1.0f/(*adt);

  for (int n=0; n<4; n++) {
//     printf("adti= %lf, res[%d] = %lf\n", adti, n, res[n]);
    del = adti*res[n];
//     printf("del is= %lf\n", del);
    q[n] = qold[n] - del;
    res[n] = 0.0f;
    *rms += del*del;
  }
}


int main(int argc, char **argv){

  long *becell, *ecell, *bound, *bedge, *edge, *cell;
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

  cell = (long *) malloc(4*ncell*sizeof(long));
  edge = (long *) malloc(2*nedge*sizeof(long));
  ecell = (long *) malloc(2*nedge*sizeof(long));
  bedge = (long *) malloc(2*nbedge*sizeof(long));
  becell = (long *) malloc( nbedge*sizeof(long));
  bound = (long *) malloc( nbedge*sizeof(long));

  x = (double *) malloc(2*nnode*sizeof(double));
  q = (double *) malloc(4*ncell*sizeof(double));
  qold = (double *) malloc(4*ncell*sizeof(double));
  res = (double *) malloc(4*ncell*sizeof(double));
  adt = (double *) malloc( ncell*sizeof(double));

  for (long n=0; n<nnode; n++) {
    if (fscanf(fp,"%lf %lf \n",&x[2*n], &x[2*n+1]) != 2) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (long n=0; n<ncell; n++) {
    if (fscanf(fp,"%ld %ld %ld %ld \n",&cell[4*n ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (long n=0; n<nedge; n++) {
    if (fscanf(fp,"%ld %ld %ld %ld \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      printf("error reading from new_grid.dat\n"); exit(-1);
    }
  }

  for (long n=0; n<nbedge; n++) {
    if (fscanf(fp,"%ld %ld %ld %ld \n",&bedge[2*n],&bedge[2*n+1],
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

  for (long n=0; n<ncell; n++) {
    for (long m=0; m<4; m++) {
        q[4*n+m] = qinf[m];
      res[4*n+m] = 0.0f;
      
    }
  }

  niter = 1000;
  printf("Before giant for loop\n");

  clock_t start, end;
  start = clock();
  
  for(int iter=1; iter<=niter; iter++) {

// save old flow solution

    for (long i = 0; i < ncell; ++i) {
      save_soln(&q[4*i], &qold[4*i]);
    }

// predictor/corrector update loop


    for(long k=0; k<2; k++) {

// calculate area/timstep

      for (long i = 0; i < ncell; ++i) {
        adt_calc(&x[2*cell[4*i]], &x[2*cell[4*i+1]], &x[2*cell[4*i+2]], &x[2*cell[4*i+3]], &q[4*i], &adt[i]);
      }
      

// calculate flux residual

      for (long i = 0; i < nedge; ++i) {
        res_calc(&x[2*edge[2*i]], &x[2*edge[2*i+1]], &q[4*ecell[2*i]], &q[4*ecell[2*i+1]], &adt[ecell[2*i]], &adt[ecell[2*i+1]], &res[4*ecell[2*i]], &res[4*ecell[2*i+1]]);
      }
      
      for (long i = 0; i < nbedge; ++i) {
        bres_calc(&x[2*bedge[2*i]], &x[2*bedge[2*i+1]], &q[4*becell[i]], &adt[becell[i]], &res[4*becell[i]], &bound[i]);
      }

// update flow field

      rms = 0.0;

      for (long i = 0; i < ncell; ++i) {
        update(&qold[4*i], &q[4*i], &res[4*i], &adt[i], &rms);
      }
    }

// print iteration history

    rms = sqrt(rms/(double) ncell);

    if (iter%100 == 0)
      printf(" %d %10.5e \n",iter,rms);
  }
  
  end = clock();
  printf("time taken: %lf seconds\n", (double)(end - start)/CLOCKS_PER_SEC);
  
}
