/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
Part of this code is taken from Numerical receipies in C, William H. Press, Saul A. Teukolsky, William T. Vetterling,
and Brian P. Flannery, CAMBRIDGE UNIVERSITY PRESS
Copyright @ Sabyasachi Sutradhar
*/


#define NR_END 1
#define FREE_ARG char*
////////////////////////////////////////

double ***dtensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
double ***t;
/* allocate pointers to pointers to rows */
t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
if (!t) printf("allocation failure 1 in f3tensor()");
t += NR_END;
t -= nrl;
/* allocate pointers to rows and set pointers to them */
t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
if (!t[nrl]) printf("allocation failure 2 in f3tensor()");
t[nrl] += NR_END;
t[nrl] -= ncl;
/* allocate rows and set pointers to them */
t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
if (!t[nrl][ncl]) printf("allocation failure 3 in f3tensor()");
t[nrl][ncl] += NR_END;
t[nrl][ncl] -= ndl;
for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
for(i=nrl+1;i<=nrh;i++) {
t[i]=t[i-1]+ncol;
t[i][ncl]=t[i-1][ncl]+ncol*ndep;
for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
}
/* return pointer to array of pointers to rows */
return t;
}
//////////////////////////////////////////////////////

int ***itensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
int ***t;
/* allocate pointers to pointers to rows */
t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
if (!t) printf("allocation failure 1 in f3tensor()");
t += NR_END;
t -= nrl;
/* allocate pointers to rows and set pointers to them */
t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
if (!t[nrl]) printf("allocation failure 2 in f3tensor()");
t[nrl] += NR_END;
t[nrl] -= ncl;
/* allocate rows and set pointers to them */
t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
if (!t[nrl][ncl]) printf("allocation failure 3 in f3tensor()");
t[nrl][ncl] += NR_END;
t[nrl][ncl] -= ndl;
for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
for(i=nrl+1;i<=nrh;i++) {
t[i]=t[i-1]+ncol;
t[i][ncl]=t[i-1][ncl]+ncol*ndep;
for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
}
/* return pointer to array of pointers to rows */
return t;
}

//////////////////////////////////////////////////////////////////////////////////
double **dvector(int nrl, int nrh, int ncl, int nch )
{
  int i;
  double **m;

  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) printf("allocation failure 1 in matrix()\n\n");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
    m[i] -= ncl;
  }
  return m;
}
//////////////////////////////////////////////

int **ivector( int nrl, int nrh, int ncl, int nch )
{
  int i;
  int **m;

  m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m) printf("allocation failure 1 in matrix()\n\n");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
    m[i] -= ncl;
  }
  return m;
}

double *dnum( int nl, int nh)
{
  double *v;

  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) printf("allocation failure in dynein()\n\n");
  return v-nl;
}
int *inum( int nl, int nh)
{
  int *v;

  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
  if (!v) printf("allocation failure in dynein()\n\n");
  return v-nl;
}
