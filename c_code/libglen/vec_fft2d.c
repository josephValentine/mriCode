#include<pthread.h>
#include<stdio.h>
#include<stdlib.h>

#define MAXTHREADS 8

/* 
  Does 2d fft on rows x cols (double) matrix with real part in re,
  imaginary part in im.  Puts result in re2 (real part) and im2
  (imaginary part).  fft2d(x,y,x,y,...) overwrites the matrix
  in x and y with its FFT.  Does inverse transform if inv != 0.
  Data is in row major format, i.e. re[0] is row 0, col 0, 
  re[1] is row 0, col 1, etc.
  If rows or cols is not a power of two, does (slow) dft on that
  dimension.
  Vectorized version with pthreads.

*/

struct threadinfo {
  double *re;
  double *im;
  int npts;
  int inv;
};
  

void *p_fft(void *stuff)
{
  struct threadinfo *s;

  s = (struct threadinfo *) stuff;
  //  printf("thread %d: starting at address (%ld,%ld) with %d points\n",pthread_self(),s->re,s->im,s->npts);
  fft(s->re,s->im,s->npts,s->inv);
  pthread_exit(NULL);
}

double *rbuf,*ibuf;
struct threadinfo ti[MAXTHREADS];
pthread_t threads[MAXTHREADS];
pthread_attr_t attr;

vec_fft2d(re,im,re2,im2,rows,cols,inv)
double *re,*im,*re2,*im2;
int rows,cols,inv;

{
  int i,j,todo,rc,status;
  //  double *blah;

  if(re != re2)
    for(i=0; i<rows*cols; i++)
      re2[i] = re[i];
  if(im != im2)
    for(i=0; i<rows*cols; i++)
      im2[i] = im[i];

/* transform rows */

/* Farm this out to threads */

/* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(i=0; i<rows;) {
    todo = (MAXTHREADS < (rows-i))?(MAXTHREADS):(rows-i);
    //    printf("vectorizing chunk of %d transforms, %d total left...\n",todo,rows - i);
    fflush(stdout);
    for(j=0; j<todo; j++,i++) {
      ti[j].re = re2+i*cols;
      ti[j].im = im2+i*cols;
      ti[j].npts = cols;
      ti[j].inv = inv;
      rc = pthread_create(&threads[j], &attr, p_fft, (void *)&ti[j]);
      if(rc) {
	fprintf(stderr,"pthread_create() failed with error code %d\n",rc);
	exit(1);
      }
    }
    for(j=0; j<todo; j++) {
      status = 1;
      rc = pthread_join(threads[j],(void **)&status);
      if(rc) {
	fprintf(stderr,"pthread_join failed for thread %d, status %s\n",j,status);
	exit(1);
      }
    }
  }

  /* transform columns */

  if((rbuf = (double *)malloc(sizeof(double)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((ibuf = (double *)malloc(sizeof(double)*rows*cols)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  //  blah = ibuf;
  for(i=0; i<rows * cols; i++) {
    rbuf[(i%cols)*rows + i/cols] = re2[i];
    ibuf[(i%cols)*rows + i/cols] = im2[i];
  }

/* Farm this out to threads */

  for(i=0; i<cols;) {
    todo = (MAXTHREADS < (cols-i))?(MAXTHREADS):(cols-i);
    //    printf("vectorizing chunk of %d transforms, %d total left...\n",todo,cols - i);
    fflush(stdout);
    for(j=0; j<todo; j++,i++) {
      ti[j].re = rbuf+i*rows;
      ti[j].im = ibuf+i*rows;
      ti[j].npts = rows;
      ti[j].inv = inv;
      //      printf("setting up thread with re = %ld, im = %ld (rbuf = %ld, ibuf = %ld, blah = %ld)\n",ti[j].re,ti[j].im,rbuf,ibuf,blah);
      rc = pthread_create(&threads[j], &attr, p_fft, (void *)&ti[j]);
      if(rc) {
	fprintf(stderr,"pthread_create() failed with error code %d\n",rc);
	exit(1);
      }
    }
    for(j=0; j<todo; j++) {
      status = 1;
      rc = pthread_join(threads[j],(void **)&status);
      if(rc) {
	fprintf(stderr,"pthread_join failed for thread %d, status %s\n",j,status);
	exit(1);
      }
    }
  }

  for(i=0; i<rows * cols; i++) {
    re2[i] = rbuf[(i%cols)*rows + i/cols];
    im2[i] = ibuf[(i%cols)*rows + i/cols];
  }
  free(ibuf);
  free(rbuf);
}

