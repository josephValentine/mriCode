#include<stdlib.h>
#include<stdio.h>
#include<errno.h>

/*
  version for floats, not doubles...
  does n point 1d fft on complex vector, real part in r, imaginary in im.
  If inv != 0, does inverse fft.
*/

f_fft(r,im,n,inv)
float *r, *im;
int n,inv;

{
  float *buf;
  int i;
  int ispow2();

  if(ispow2(n)) {
    if((buf = (float *)malloc(sizeof(float) * n * 2)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    for(i=0; i<n; i++) {
      buf[2*i] = r[i];
      buf[2*i+1] = im[i];
    }
    f_four1(buf-1,n,(inv)?(-1):(1));
    if(inv) 
      for(i=0; i<n; i++) {
	r[i] = buf[2*i] / (float)n;
	im[i] = buf[2*i+1] / (float)n;
      }
    else
      for(i=0; i<n; i++) {
	r[i] = buf[2*i];
	im[i] = buf[2*i+1];
      }
    free(buf);
  }
  else {
    f_dft(r,im,n,inv);
  }
}

