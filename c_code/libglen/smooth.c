/* optimally fits function contained in double pointers "inr" and "ini" 
   (real and imaginary)
   of size "ydim" by "xdim" over set of points defined by costfunction
   "cost" (doubles, also ydim by xdim)
   using "yres" by "xres" fourier coefficients.  
   "yres" and "xres" are assumed to be less than "ydim" and "xdim".
   If outr and outi are null
   pointers, result is put in inr and ini.  If ini is null pointer, fit is
   assumed to be purely real.
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

smooth(inr,ini,ydim,xdim,outr,outi,yres,xres,cost)
double *inr,*ini,*outr,*outi;
int ydim,xdim,yres,xres;
double *cost;

{
  double *Ar,*Ai,*Fr,*Fi,**Sr,**Si,*br,*bi,ddum,max;
  int i,j,ui,vi,uj,vj,fold,ncoeff,*index,real;
  FILE *debug;

  debug = fopen("debugout","w");

  real = 0;

/* resolution should always be odd */

  xres -= xres % 2;
  xres++;
  yres -= yres % 2;
  yres++;

  ncoeff = xres * yres;
  printf("ncoeff = %d\n",ncoeff);
   if((Sr = (double **)malloc(sizeof(double *) *ncoeff)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  for(i=0; i<ncoeff; i++)
    if((Sr[i] = (double *)malloc(sizeof(double) * ncoeff)) == NULL) {
      perror("malloc failed...");
      exit(1);
    }
  if((Si = (double **)malloc(sizeof(double *) *ncoeff)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  for(i=0; i<ncoeff; i++)
    if((Si[i] = (double *)malloc(sizeof(double) * ncoeff)) == NULL) {
      perror("malloc failed...");
      exit(1);
    }
  if((br = (double *)malloc(sizeof(double) * ncoeff)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((bi = (double *)malloc(sizeof(double) * ncoeff)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((index = (int *)malloc(sizeof(int) * ncoeff)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((Ar = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((Ai = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((Fr = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  if((Fi = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
    perror("malloc failed...");
    exit(1);
  }
  fold = (outr == (double *)0);
  if(fold)
    if((outr = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
      perror("malloc failed...");
      exit(1);
    }
  if(outi == (double *)0) {
    real = 1;
    if((outi = (double *)malloc(sizeof(double) * xdim * ydim)) == NULL) {
      perror("malloc failed...");
      exit(1);
    }
  }

  for(i=0; i<xdim*ydim; i++) {
    Ar[i] = cost[i];
    Ai[i] = 0.0;
    Fr[i] = inr[i] * cost[i];
    if(ini)
      Fi[i] = ini[i] * cost[i];
    else
      Fi[i] = 0.0;
  }
  printf("transforming f...");
  fflush(stdout);
  fft2d(Fr,Fi,Fr,Fi,ydim,xdim,0);
  printf("done\n");
  printf("transforming a...");
  fflush(stdout);
  fft2d(Ar,Ai,Ar,Ai,ydim,xdim,0);
  printf("done\n");
  for(i=0; i<ncoeff; i++) {
    for(j=0; j<ncoeff; j++) {
      ui = -xres/2 + i / yres;
      vi = -yres/2 + i % yres;
      uj = -xres/2 + j / yres;
      vj = -yres/2 + j % yres;
      Sr[j][i] = Ar[((uj-ui+xdim)%xdim) + ((vj-vi+ydim)%ydim) *xdim];
      Si[j][i] = Ai[((uj-ui+xdim)%xdim) + ((vj-vi+ydim)%ydim) *xdim];
    }
  }
  free(Ar);
  free(Ai);
  for(i=0; i<ncoeff; i++) {
    ui = -xres/2 + i / yres;
    vi = -yres/2 + i % yres;
    br[i] = Fr[(ui+xdim)%xdim + ((vi + ydim)%ydim)*xdim];
    bi[i] = Fi[(ui+xdim)%xdim + ((vi + ydim)%ydim)*xdim];
  }
  
  printf("done\n");
  
  for(i=0; i<ncoeff; i++) {  /* lu decomp routine uses vectors starting with x[1] */
    Sr[i]--;
    Si[i]--;
  }
  printf("decomposing...");
  fflush(stdout);
  cludcmp(Sr-1,Si-1,ncoeff,index-1,&ddum);
  printf("done\nsolving...");
  fflush(stdout);
  clubksb(Sr-1,Si-1,ncoeff,index-1,br-1,bi-1);
  printf("done\n");
  free(Sr);
  free(Si);

  for(i=0; i<xdim; i++)
    for(j=0; j<ydim; j++) {
      Fr[i + j * xdim] = 0.0;
      Fi[i + j * xdim] = 0.0;
    }
  for(i=0; i<ncoeff; i++) {
    ui = -xres/2 + i / yres;
    vi = -yres/2 + i % yres;
    ui = -ui;  /* don't ask me why, just tweaking, I'll figure it out later. */
    vi = -vi;
    Fr[(ui+xdim)%xdim + ((vi + ydim)%ydim)*xdim] = br[i];
    Fi[(ui+xdim)%xdim + ((vi + ydim)%ydim)*xdim] = bi[i];
  }
  fwrite(Fr,sizeof(double),xdim*ydim,debug);
  fwrite(Fi,sizeof(double),xdim*ydim,debug);
  free(br);
  free(bi);
  printf("inverse transforming...");
  fflush(stdout);
  if(fold)
    fft2d(Fr,Fi,inr,ini,ydim,xdim,0);
  else
    fft2d(Fr,Fi,outr,outi,ydim,xdim,0);
  if(real) {
    max = 0.0;
    for(i=0; i<xdim*ydim; i++)
      if(fabs(outi[i]) > max)
	max = fabs(outi[i]);
    printf("smooth: max. imag. value = %g\n",max);
    fwrite(outi,sizeof(double),xdim*ydim,debug);
  }
  free(Fr);
  free(Fi);
  printf("done\n");
}
  
  






