#include<stdlib.h>
#include<stdio.h>

/* 
  Version for floats, not doubles...

  Does 2d fft on rows x cols (float) matrix with real part in re,
  imaginary part in im.  Puts result in re2 (real part) and im2
  (imaginary part).  fft2d(x,y,x,y,...) overwrites the matrix
  in x and y with its FFT.  Does inverse transform if inv != 0.
  Data is in row major format, i.e. re[0] is row 0, col 0, 
  re[1] is row 0, col 1, etc.
*/

f_fft2d(re,im,re2,im2,rows,cols,inv)
float *re,*im,*re2,*im2;
int rows,cols,inv;

{
  int i;
  float *rbuf,*ibuf;

  if(re != re2)
    for(i=0; i<rows*cols; i++)
      re2[i] = re[i];
  if(im != im2)
    for(i=0; i<rows*cols; i++)
      im2[i] = im[i];

/* transform rows */

  for(i=0; i<rows; i++)
    f_fft(re2+i*cols,im2+i*cols,cols,inv);
  
  /* transform columns */

  if((rbuf = (float *)malloc(sizeof(float)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((ibuf = (float *)malloc(sizeof(float)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i=0; i<rows * cols; i++) {
    rbuf[(i%cols)*rows + i/cols] = re2[i];
    ibuf[(i%cols)*rows + i/cols] = im2[i];
  }
  for(i=0; i<cols; i++)
    f_fft(rbuf+i*rows,ibuf+i*rows,rows,inv);
  for(i=0; i<rows * cols; i++) {
    re2[i] = rbuf[(i%cols)*rows + i/cols];
    im2[i] = ibuf[(i%cols)*rows + i/cols];
  }
  free(ibuf);
  free(rbuf);
}

