#include<stdlib.h>
#include<stdio.h>

/* 
  Does 3d fft on rows x cols x planes (double) matrix with real part in re,
  imaginary part in im.  Puts result in re2 (real part) and im2
  (imaginary part).  fft3d(x,y,x,y,...) overwrites the matrix
  in x and y with its FFT.  Does inverse transform if inv != 0.
  Data is in row major format, i.e. re[0] is row 0, col 0, plane 0
  re[1] is row 0, col 1, plane 0, re[row*col] is row 0, col 0, plane1, etc.
*/




fft3d(re,im,re2,im2,rows,cols,planes,inv)
double *re,*im,*re2,*im2;
int rows,cols,planes,inv;

{
  int i,j,k;
  double *rbuf,*ibuf;
  int a;

  if(re != re2)
    for(i=0; i<rows*cols*planes; i++)
      re2[i] = re[i];
  if(im != im2)
    for(i=0; i<rows*cols*planes; i++)
      im2[i] = im[i];

/* transform planes */

  a = cols * rows;
  for(i=0; i<planes; i++)
    fft2d(re2+i*a,im2+i*a,re2+i*a,im2+i*a,rows,cols,inv);
  
/* transform in 3rd dimension */

  if((rbuf = (double *)malloc(sizeof(double)*planes)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((ibuf = (double *)malloc(sizeof(double)*planes)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0; i<rows; i++)
    for(j=0; j<cols; j++) {
      for(k=0; k<planes; k++) {
	rbuf[k] = re2[k*a+i*cols+j];
	ibuf[k] = im2[k*a+i*cols+j];
      }
      fft(rbuf,ibuf,planes,inv);
      for(k=0; k<planes; k++) {
	re2[k*a+i*cols+j] = rbuf[k];
	im2[k*a+i*cols+j] = ibuf[k];
      }
    }
  free(ibuf);
  free(rbuf);
}
