#include<stdlib.h>
#include<stdio.h>
#include<errno.h>

/*
  does n point 1d fft on complex vector, real part in r, imaginary in im.
  If inv != 0, does inverse fft.  If n is not a power of two, does (slow) dft.
*/

fft(r,im,n,inv)
double *r, *im;
int n,inv;

{
  double *buf;
  int i;
  int ispow2();

  if(ispow2(n)) {
    if((buf = (double *)malloc(sizeof(double) * n * 2)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    for(i=0; i<n; i++) {
      buf[2*i] = r[i];
      buf[2*i+1] = im[i];
    }
    four1(buf-1,n,(inv)?(-1):(1));
    if(inv) 
      for(i=0; i<n; i++) {
	r[i] = buf[2*i] / (double)n;
	im[i] = buf[2*i+1] / (double)n;
      }
    else
      for(i=0; i<n; i++) {
	r[i] = buf[2*i];
	im[i] = buf[2*i+1];
      }
    free(buf);
  }
  else {
    dft(r,im,n,inv);
  }
}

ispow2(n)
int n;

{
  if(n < 2)
    return(0);
  for(;!(n%2);n /= 2);
  return(((n == 1)?(1):(0)));
}
