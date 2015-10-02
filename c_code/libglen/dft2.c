/*
  Discrete fourier transform routine for doubles.
  Not fast fourier transform.  Works for any n, not just powers of two,
  but takes a long time.  Inverse transform calculated if inv != 0.
  Note that fft.c calls dft if n is not power of two, so just use routine
  fft() for everything.
  Trying static pointers to sin and cosin tables to make faster for repeated calls with same vector length
  This is SORT OF thread-safe, in that only one thread can initialize the size of the sin table at a time,
  but if a later thread decides to change it, threads already in the function will suffer.  Presumably,
  we are dealing with only one "n" for all threads.  If so, this is thread-safe.
*/

#include<stdlib.h>
#include<math.h>
#include<stdio.h>

//extern pthread_mutex_t mutex;

dft2(re,im,n,inv)
double *re,*im;
int n,inv;

{
  double sdum,cdum;
  double *fr,*fi;
  int i,j,k;
  static int lasttime = 0;
  static double *sintab = 0;
  static double *costab = 0;

  if(n != lasttime) {
    //    pthread_mutex_lock(&mutex);
    if(n != lasttime) {  // check again, could have changed while waiting for mutex
      lasttime = n;
      if(sintab)
	free(sintab);
      if(costab)
	free(costab);
      if((sintab = (double *)malloc(sizeof(double)*n)) == NULL) {
	fprintf(stderr,"dft: malloc failed.\n");
	return(-1);
      }
      if((costab = (double *)malloc(sizeof(double)*n)) == NULL) {
	fprintf(stderr,"dft: malloc failed.\n");
	free(sintab);
	return(-1);
      }
      for(i=0; i<n; i++) {
	sintab[i] = sin(2.0 * M_PI * (double)i / (double)n);
	costab[i] = cos(2.0 * M_PI * (double)i / (double)n);
      }
    }
    //    pthread_mutex_unlock(&mutex);
  }
  if((fr = (double *)malloc(sizeof(double)*n)) == NULL) {
    fprintf(stderr,"dft: malloc failed.\n");
    return(-1);
  }
  if((fi = (double *)malloc(sizeof(double)*n)) == NULL) {
    fprintf(stderr,"dft: malloc failed.\n");
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
    for(j=1,k=i; j<n; j++) {
      if(inv) {
	re[i] += fr[j] * costab[k] + fi[j] * sintab[k];
	im[i] += -fr[j] * sintab[k] + fi[j] * costab[k];
      }
      else {
	re[i] += fr[j] * costab[k] - fi[j] * sintab[k];
	im[i] += fr[j] * sintab[k] + fi[j] * costab[k];
      }
      k = (k + i)%n;
    }
  }
  if(inv)
    for(i=0;i<n;i++) {
      re[i] /= (double)n;
      im[i] /= (double)n;
    }
  free(fr);
  free(fi);
  return(0);
}

      
