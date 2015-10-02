/*
  Discrete fourier transform routine for floats.
  Not fast fourier transform.  Works for any n, not just powers of two,
  but takes a long time.  Inverse transform calculated if inv != 0.
  Since f_fft calls f_dft if n is not power of two, just use f_fft for 
  everything.
*/

#include<stdlib.h>
#include<math.h>
#include<stdio.h>

f_dft(re,im,n,inv)
float *re,*im;
int n,inv;

{
  double theta,thinc;
  double sdum,cdum;
  float *fr,*fi;
  int i,j;

  if((fr = (float *)malloc(sizeof(float)*n)) == NULL) {
    fprintf(stderr,"f_dft: malloc failed.\n");
    return(-1);
  }
  if((fi = (float *)malloc(sizeof(float)*n)) == NULL) {
    fprintf(stderr,"f_dft: malloc failed.\n");
    free(fr);
    return(-1);
  }
  for(i=0;i<n;i++) {
    fr[i] = re[i];
    fi[i] = im[i];
  }
  re[0] = im[0] = 0.0;
  for(j=0;j<n;j++) {
    re[0] += fr[j];
    im[0] += fi[j];
  }
  for(i=1;i<n;i++) {
    re[i] = fr[0];
    im[i] = fi[0];
    thinc = ((inv)?(-1.0):(1.0)) * 2.0 * M_PI * (double)i / (double)n;
    theta = thinc;
    for(j=1; j<n; j++) {
      sdum = sin(theta);
      cdum = cos(theta);
      re[i] += fr[j] * cdum - fi[j] * sdum;
      im[i] += fr[j] * sdum + fi[j] * cdum;
      theta += thinc;
    }
  }
  if(inv)
    for(i=0;i<n;i++) {
      re[i] /= (float)n;
      im[i] /= (float)n;
    }
  free(fr);
  free(fi);
  return(0);
}

      
