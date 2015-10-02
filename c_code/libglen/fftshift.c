#include<stdlib.h>
#include<stdio.h>

/* 
  shift origin of vector of size n, with n even
*/

fftshift(vec,n)
     double *vec;
     int n;
{
  double *dum;
  int i,org;

  if((dum = (double *)malloc(sizeof(double)*n)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0; i<n; i++)
    dum[i] = vec[i];
  org = n/2;
  for(i=0; i<n; i++)
    vec[i] = dum[(i + org)%n];
  free(dum);
}

/*
  shift origin of matrix
*/

fftshift2d(mat,rows,cols)
double *mat;
int rows,cols;

{
  double *dum;
  int i,j,xorg,yorg,pts;

  pts = rows*cols;
  if((dum = (double *)malloc(sizeof(double)*pts)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0;i<pts;i++)
    dum[i] = mat[i];
  xorg = cols/2;
  yorg = rows/2;
  for(i=0;i<cols;i++)
    for(j=0;j<rows;j++) 
      mat[i + j*cols] = dum[(i+xorg)%cols+((j + yorg)%rows)*cols];
  free(dum);
}

/*
  like above, but shifts one more pixel if array sizes are odd.  
  use to shift from upper left corner to center.  To go from center
  to upper left corner, use fftshift2d.
*/
  
oddfftshift2d(mat,rows,cols)
double *mat;
int rows,cols;

{
  double *dum;
  int i,j,xorg,yorg,pts;

  pts = rows*cols;
  if((dum = (double *)malloc(sizeof(double)*pts)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0;i<pts;i++)
    dum[i] = mat[i];
  if((cols % 2))
    xorg = cols/2 + 1;
  else
    xorg = cols/2;
  if((rows % 2))
    yorg = rows/2 + 1;
  else
    yorg = rows/2;
  for(i=0;i<cols;i++)
    for(j=0;j<rows;j++) 
      mat[i + j*cols] = dum[(i+xorg)%cols+((j + yorg)%rows)*cols];
  free(dum);
}

/* float versions of above routines .. */

f_fftshift2d(mat,rows,cols)
float *mat;
int rows,cols;

{
  float *dum;
  int i,j,xorg,yorg,pts;

  pts = rows*cols;
  if((dum = (float *)malloc(sizeof(float)*pts)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0;i<pts;i++)
    dum[i] = mat[i];
  xorg = cols/2;
  yorg = rows/2;
  for(i=0;i<cols;i++)
    for(j=0;j<rows;j++) 
      mat[i + j*cols] = dum[(i+xorg)%cols+((j + yorg)%rows)*cols];
  free(dum);
}

f_oddfftshift2d(mat,rows,cols)
float *mat;
int rows,cols;

{
  float *dum;
  int i,j,xorg,yorg,pts;

  pts = rows*cols;
  if((dum = (float *)malloc(sizeof(float)*pts)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0;i<pts;i++)
    dum[i] = mat[i];
  if((cols % 2))
    xorg = cols/2 + 1;
  else
    xorg = cols/2;
  if((rows % 2))
    yorg = rows/2 + 1;
  else
    yorg = rows/2;
  for(i=0;i<cols;i++)
    for(j=0;j<rows;j++) 
      mat[i + j*cols] = dum[(i+xorg)%cols+((j + yorg)%rows)*cols];
  free(dum);
}

  

