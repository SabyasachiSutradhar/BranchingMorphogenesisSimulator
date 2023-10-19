/*
This header file written by Sabyasachi Sutradhar, Yale university (sabyasachi.Sutradhar@yale.edu) 2022
This code simulates branching morphogensis of Drosophila class-IV dendritic arbor. branches are generates from the origin (0,0)
and then the tips go through three state dynamics G-P-S with different rates and also new branches spawn from an existing branchesrandomly
with a rate.
Part of this code is taken from Numerical receipies in C, William H. Press, Saul A. Teukolsky, William T. Vetterling,
and Brian P. Flannery, CAMBRIDGE UNIVERSITY PRESS
Copyright @ Sabyasachi Sutradhar
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NR_END 1
#define FREE_ARG char*

double sqr(double q){
  return (q*q);
}
double cube(double p){
  return(p*p*p);
}
double cube_root(double p){
  return(pow(p,1.0/3.0));
}
int find_min(int i,int j){
  if(i>=j) return(j);
  else return(i);
}
int find_max(int i,int j){
  if(i>=j) return(i);
  else return(j);
}
int find_dmin(double i,double j){
  if(i<j) return(i);
  else return(j);
}
int find_dmax(double i,double j){
  if(i>=j) return(i);
  else return(j);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(char const error_text[])
/* Numerical Recipes standard error handler */ {
printf("Numerical Recipes run-time error...\n");
printf("%s\n",error_text);
printf("...now exiting to system...\n");
exit(1);
}
float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */ {
float *v;
v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float))); if (!v) nrerror("allocation failure in vector()");
return v-nl+NR_END;
}
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */ {
int *v;
v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))); if (!v) nrerror("allocation failure in ivector()"); return v-nl+NR_END;
}
unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */ {
unsigned char *v;
v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char))); if (!v) nrerror("allocation failure in cvector()");
return v-nl+NR_END;
}
unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */ {
unsigned long *v;
v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long))); if (!v) nrerror("allocation failure in lvector()");
return v-nl+NR_END;
}
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */ {
double *v;
v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); if (!v) nrerror("allocation failure in dvector()");
return v-nl+NR_END;
}
float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; float **m;
/* allocate pointers to rows */
m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
if (!m) nrerror("allocation failure 1 in matrix()"); m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float))); if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
 return m; }
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */ {
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; double **m;
/* allocate pointers to rows */
m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*))); if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double))); if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m; }
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */ {
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1; int **m;
/* allocate pointers to rows */
m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*))); if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int))); if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m; }
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */ {
long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl; float **m;
/* allocate array of pointers to rows */
m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*))); if (!m) nrerror("allocation failure in submatrix()");
m += NR_END;
m -= newrl;
/* set pointers to rows */ for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
/* return pointer to array of pointers to rows */
return m; }
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch) /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1 and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1; float **m;
/* allocate pointers to rows */
m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*))); if (!m) nrerror("allocation failure in convert_matrix()"); m += NR_END;
m -= nrl;
/* set pointers to rows */
m[nrl]=a-ncl;
for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol; /* return pointer to array of pointers to rows */ return m;
}
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) /* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1; float ***t;
/* allocate pointers to pointers to rows */
t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**))); if (!t) nrerror("allocation failure 1 in f3tensor()");
t += NR_END;
t -= nrl;
/* allocate pointers to rows and set pointers to them */
t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*))); if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
t[nrl] += NR_END;
t[nrl] -= ncl;
/* allocate rows and set pointers to them */
t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float))); if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
t[nrl][ncl] += NR_END; t[nrl][ncl] -= ndl;
for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep; for(i=nrl+1;i<=nrh;i++) {
t[i]=t[i-1]+ncol; t[i][ncl]=t[i-1][ncl]+ncol*ndep; for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
}
/* return pointer to array of pointers to rows */
return t; }
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */ {
free((FREE_ARG) (v+nl-NR_END)); }
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) /* free a float matrix allocated by matrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END)); }
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) /* free a double matrix allocated by dmatrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END)); }
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) /* free an int matrix allocated by imatrix() */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END)); }
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch) /* free a submatrix allocated by submatrix() */
{
free((FREE_ARG) (b+nrl-NR_END)); }
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch) /* free a matrix allocated by convert_matrix() */
{
free((FREE_ARG) (b+nrl-NR_END)); }
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */ {
free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END)); free((FREE_ARG) (t[nrl]+ncl-NR_END)); free((FREE_ARG) (t+nrl-NR_END));
}
double distance2D(double xi,double yi,double xf,double yf){
 return(sqrt((xi-xf)*(xi-xf)+(yi-yf)*(yi-yf)));
}
double distance(double xi,double yi,double zi ,double xf,double yf,double zf){
 return(sqrt(sqr(xi-xf)+sqr(yi-yf)+sqr(zi-zf)));
}
double random_sign(){
  if(ran2(&idum)>0.5) return(1.0);
  else return(-1.0);
}
double Volume(double r){
  return 4.0*pi*r*r*r/3.0;
}
double Surface_area(double r){
  return 4.0*pi*r*r;
}
double x_within_circle(double a){
  if(ran2(&idum)>0.5) return(ran2(&idum)*a);
  else return (-ran2(&idum)*a);
}
double y_within_circle(double x,double a,double b){
  if(ran2(&idum)>0.5)   return(ran2(&idum)*b*sqrt(1-sqr(x)/sqr(a)));
  else  return(-ran2(&idum)*b*sqrt(1-sqr(x)/sqr(a)));
}
double z_within_circle(double x,double a,double y,double b,double c){
  if(ran2(&idum)>0.5)  return(ran2(&idum)*c*sqrt(1-(sqr(x)/sqr(a)+sqr(y)/sqr(b))));
  else  return(-ran2(&idum)*c*sqrt(1-(sqr(x)/sqr(a)+sqr(y)/sqr(b))));
}
double wrap_angle(double theta){
  if(theta<0.0) {theta=2.0*pi+theta;}
  if(theta>2.0*pi){theta=fmod(theta,2.0*pi);}
  return(theta);
}
int random_int(int a,int b){
  return (a+floor(ran2(&idum)*(b-a+1)));
}

int lineseg_intersection(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
int ind=0;
double den=(x4-x3)*(y1-y2)-(x1-x2)*(y4-y3);
if(den!=0.0000){
double ta,tb;
ta=((y3-y4)*(x1-x3)+(x4-x3)*(y1-y3))/den;
tb=((y1-y2)*(x1-x3)+(x2-x1)*(y1-y3))/den;
if(ta>=0 && ta<=1){
if(tb>=0 && tb<=1){
ind=1;
}
}
}
return ind;
}

double **equispace_points(double *x,double * y,int N1,int N2){
  double **p,**q;
  p=dmatrix(0,N1,0,2);
  q=dmatrix(0,N2,0,2);
  for(int i=1;i<=N1;i++){
    p[i][1]=x[i];
    p[i][2]=y[i];
  }
  double *currentpt,totaldist,intv,disttmp,distsum,*ptnow,*pttarget,remainder,*newpt;
  currentpt=dvector(0,2);
  ptnow=dvector(0,2);
  pttarget=dvector(0,2);
  newpt=dvector(0,2);

  double l,m,norm;
  int indfirst,len,kk,ind;
  indfirst=2;

  for(unsigned int j=1;j<=2;j++){
    currentpt[j]=p[1][j];
    q[1][j]=currentpt[j];
  }

  totaldist=0.0;
  for(int i=1;i<=N1-1;i++){
    totaldist+=distance2D(p[i][1],p[i][2],p[i+1][1],p[i+1][2]);
  }
  intv=totaldist/(double)(N2-1);
///////////////////////////////////////////////
  for(int k=1;k<=N2-1;k++){
//printf("%d %d\n",k,indfirst);
    for(int j=1;j<=2;j++){
      ptnow[j]=currentpt[j];
      pttarget[j]= p[indfirst][j];
    }
    remainder = intv;
    kk = 0;
    ind=0;
    distsum = 0.0;
    while(ind==0){
      disttmp=distance2D(ptnow[1],ptnow[2],pttarget[1],pttarget[2]);
      distsum  = distsum+disttmp;
      if(distsum >= intv){
        norm=distance2D(ptnow[1],ptnow[2],pttarget[1],pttarget[2]);
        l=(pttarget[1]-ptnow[1])/norm;
        m=(pttarget[2]-ptnow[2])/norm;
        newpt[1]=remainder*l+ptnow[1];newpt[2]=remainder*m+ptnow[2];
        ind=1;
      }else{
        remainder = remainder - disttmp;
        for(int j=1;j<=2;j++){ptnow[j] = pttarget[j];}
        kk=kk+1;
        if(indfirst+kk>N1){
          for(int j=1;j<=2;j++){newpt[j] = p[N1][j];}
                  ind=1;
        }else{
          for(int j=1;j<=2;j++){pttarget[j] = p[indfirst+kk][j];}
        }
      }
    }

    for(int j=1;j<=2;j++){
      q[k+1][j]=newpt[j];
      currentpt[j] = newpt[j];
    }
  indfirst = indfirst + kk;
  }
return q;
free_dmatrix(p,0,N1,0,2);
free_dmatrix(q,0,N2,0,2);
free_dvector(currentpt,0,2);
free_dvector(ptnow,0,2);
free_dvector(pttarget,0,2);
free_dvector(newpt,0,2);
}
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(double px, double py,double qx,double qy,double rx,double ry)
{
  double val = (qy - py) * (rx - qx) -
            (qx - px) * (ry - qy);

  if (val == 0) return 0;  // colinear
  return (val > 0)? 1: 2; // clock or counterclock wise
}

double ConvexHullArea(double *x,double *y, int n){
  double **hull;
  hull=dmatrix(0,n,0,2);

int n_hull=0;
  // Find the leftmost point
  int l = 1;
  for (int i = 2; i <= n; i++){
    if (x[i] < x[l]){l = i;}
  }
  // Start from leftmost point, keep moving counterclockwise
  // until reach the start point again.  This loop runs O(h)
  // times where h is number of points in result or output.
  int p, q;
p=l;
  do
  {
    // Add current point to result
    n_hull++;
    hull[n_hull][1]=x[p];hull[n_hull][2]=y[p];

    q =(p+1);
    if(q>n){q=q-n;}

    for (int i = 1; i <= n; i++){
      if (orientation(x[p],y[p],x[i],y[i],x[q],y[q]) == 2){
        q = i;
      }
    }

    p = q;

  } while (p != l);  // While we don't come to first point
  // calculate hull area
  double area=0.0;
  int j;
  for(int i=1;i<=n_hull;i++){
    j=i+1;
    if(j>n_hull){j=1;}
    area+=0.5*fabs(hull[i][1]*hull[j][2]-hull[j][1]*hull[i][2]);
  }
free_dmatrix(hull,0,n,0,2);
return area;
}
