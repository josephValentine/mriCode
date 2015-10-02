#include<stdlib.h>
#include<stdio.h>

/* 
  Does 2d fft on rows x cols (double) matrix with real part in re,
  imaginary part in im.  Puts result in re2 (real part) and im2
  (imaginary part).  fft2d(x,y,x,y,...) overwrites the matrix
  in x and y with its FFT.  Does inverse transform if inv != 0.
  Data is in row major format, i.e. re[0] is row 0, col 0, 
  re[1] is row 0, col 1, etc.
  If rows or cols is not a power of two, does (slow) dft on that
  dimension.
*/

fft2d(re,im,re2,im2,rows,cols,inv)
double *re,*im,*re2,*im2;
int rows,cols,inv;

{
  int i;
  double *rbuf,*ibuf;

  if(re != re2)
    for(i=0; i<rows*cols; i++)
      re2[i] = re[i];
  if(im != im2)
    for(i=0; i<rows*cols; i++)
      im2[i] = im[i];

  if((rows == 1) || (cols == 1)) { // just a 1D transform
    fft(re2,im2,(cols > rows)?(cols):(rows),inv);
    return;
  }

/* transform rows */

  for(i=0; i<rows; i++)
    fft(re2+i*cols,im2+i*cols,cols,inv);

  /* transform columns */

  if((rbuf = (double *)malloc(sizeof(double)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((ibuf = (double *)malloc(sizeof(double)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0; i<rows * cols; i++) {
    rbuf[(i%cols)*rows + i/cols] = re2[i];
    ibuf[(i%cols)*rows + i/cols] = im2[i];
  }
  for(i=0; i<cols; i++)
    fft(rbuf+i*rows,ibuf+i*rows,rows,inv);
  for(i=0; i<rows * cols; i++) {
    re2[i] = rbuf[(i%cols)*rows + i/cols];
    im2[i] = ibuf[(i%cols)*rows + i/cols];
  }
  free(ibuf);
  free(rbuf);
}

