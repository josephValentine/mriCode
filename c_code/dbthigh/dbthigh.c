// 3 point dixon recon.  Takes directory name or names as argument, looks for specified series);
// in DICOM files, tries to find three echoes to do a recon with.
// usage: 3pt series dir1 [dir2 ...]
//      series is comma separated list of series, e.g. "1,5,7" (no spaces)
// Version to do fat volume calculation for Hopkins metabolic syndrome grant

#include<stdlib.h>
#include<stdio.h>
#include<errno.h>
#include<limits.h>
#include<math.h>
#include"extern.h"
#include<complex.h>
#include<string.h>

#define SLOP 0.01 // possible fp error in comparing echo time differences to be ignored.
#define min(a,b) (((a)<(b))?(a):(b))

int ispow2();  // this function is in libglen.a

//FILE *debugfp;


main(argc,argv)
     int argc;
     char **argv;
{
  struct bean *bp;
  int i;
  FILE *fp;
  double meanfwfrac,totalpts;
  char fname[1024];

  //  debugfp = fopen("./debug","w");
  init(argc,argv);
  magphase();
  for(firstb0done=0,i=nbeans/2; i < nbeans; i++)
    doitman(i);
  for(i=(nbeans/2)-1; i>=0; i--)
    doitman(i);
  sprintf(fname,"%s/%s\0",resultsdir,"info.txt");
  fp = fopen(fname,"a");
  if(basename) {
    if(!strcmp(basename,"abd"))
      fprintf(fp,"\nAbdomen:\n");
  }
  else
    fprintf(fp,"\nLegs:\n");
  fprintf(fp,"\npatient name: %s\n",thispatient);
  fprintf(fp,"image size %d rows, %d columns\n",beanstar->rows,beanstar->cols);
  fprintf(fp,"pixel dimensions: %g mm by %g mm\n",pixelx,pixely);
  fprintf(fp,"\t(pixel area = %g mm^2, %g cm^2)\n",pixelx*pixely,pixelx*pixely/100.0);
  if(poff)
    fprintf(fp,"phase offset %g degrees\n",poff);
  totalpts = 0.0;
  meanfwfrac = 0.0;
  for(i=1; i<=nbeans; i++) {
    for(bp = beanstar; bp; bp = bp->next)
      if(bp->order == i)
	break;
    if(!bp) {
      fprintf(stderr,"ran out of beans! yikes!\n");
      exit(1);
    }
    //    fprintf(fp,"slice %d, location %g, %d valid points, %g%% average f/w frac\n",i,bp->sliceloc,bp->nvalidpts,bp->meanfwfrac/10.0);
    meanfwfrac += (bp->nvalidpts * bp->meanfwfrac);
    totalpts += (double)bp->nvalidpts;
  }
  //  fprintf(fp,"\nTotal valid points: %g  Overall fat/water fraction: %g%%\n",totalpts,meanfwfrac / totalpts / 10.0);
  fclose(fp);
}


magphase()
{
  struct image *ip1,*ip2;
  int i,done;
  struct bean *bp,*newbean;
  double sloc;

  // combine phase and magnitude images

  for(;;) {
    if(basep == (struct image *)0)  // could happen if no matches for any image
      break;
    for(ip1 = basep; ip1->done && ip1->next; ip1 = ip1->next);

    if(ip1->done)
      break;
    if(ip1->imtype != 'P')  // ignore images not of type 'P' or 'M'
      if(ip1->imtype != 'M') {
	deleteim(ip1);
	continue;
      }
    done = 0;
    for(ip2 = ip1->next;ip2;ip2 = ip2->next) {
      if((ip1->rows != ip2->rows) || (ip1->cols != ip2->cols) 
	 || strcmp(ip1->sliceloc,ip2->sliceloc) 
	 || strcmp(ip1->imagepos,ip2->imagepos)|| (ip1->accession != ip2->accession)
	 || strcmp(ip1->pixelspacing,ip2->pixelspacing)
	 || strcmp(ip1->slicethickness,ip2->slicethickness) || (ip1->tr != ip2->tr)
	 || strcmp(ip1->acqtime,ip2->acqtime) || (ip1->te != ip2->te))
	continue;
      switch(ip1->imtype) {
      case 'M':
	if(ip2->imtype != 'P')
	  continue;
	done++;
	break;
      case 'P':
	if(ip2->imtype != 'M')
	  continue;
	done++;
	break;
      default:
	fprintf(stderr,"unknown image type \"%c\" from file %s, bailing...\n",ip1->imtype,ip1->filename);
	exit(1);
      }
      if(done)
	break;
    }
    if((ip2 == (struct image *)0) && !done) {
      if(ip1->imtype == 'M')
	printf("No phase image for magnitude image %s, series %d image %d, skipping...\n",ip1->filename,ip1->series,ip1->imageno);
      else
	printf("No magnitude image for phase image %s, series %d image %d, skipping...\n",ip1->filename,ip1->series,ip1->imageno);
      deleteim(ip1);
      continue;
    }
    //    printf("matching series %d image %d to series %d image %d\n",ip1->series,ip1->imageno,ip2->series,ip2->imageno);
    if((ip1->real = (double *)malloc(sizeof(double)*ip1->rows*ip1->cols)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    if((ip1->imag = (double *)malloc(sizeof(double)*ip1->rows*ip1->cols)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    for(i=0; i<ip1->rows * ip1->cols; i++)
      switch(ip1->imtype) {
      case 'M':
	// Siemens phase images represent 0-2pi with range of 0-4096.
	ip1->real[i] = (double)ip1->image[i] * cos((double)ip2->image[i] / 4096.0 * 2.0 * M_PI);
	ip1->imag[i] = (double)ip1->image[i] * sin((double)ip2->image[i] / 4096.0 * 2.0 * M_PI);
	break;
      case 'P':
	ip1->real[i] = (double)ip2->image[i] * cos((double)ip1->image[i] / 4096.0 * 2.0 * M_PI);
	ip1->imag[i] = (double)ip2->image[i] * sin((double)ip1->image[i] / 4096.0 * 2.0 * M_PI);
	break;
      default:
	fprintf(stderr,"The sky is falling!\n");
	exit(1);
      }
    ip1->done++;
    deleteim(ip2);
    free(ip1->image);
    ip1->image = (short *)0;
  }
  if(basep == (struct image *)0) {
    printf("No images.\n");
    exit(0);
  }
  for(ip1 = basep; ip1; ip1 = ip1->next) {
    //    printf("Complex image at location %s\n",ip1->sliceloc);
    if(!byimagenumber)  // if ordering by image number, sliceloc becomes image number.
      sloc = atof(ip1->sliceloc);
    else
      sloc = ip1->imageno;
    if(!beanstar) {
      if((beanstar = (struct bean *)malloc(sizeof(struct bean))) == NULL) {
	perror("malloc failed");
	exit(1);
      }
      beanstar->sliceloc = sloc;
      beanstar->next = beanstar->prev = (struct bean *)0;
    }
    else {
      for(bp = beanstar; bp; bp = bp->next) {
	if(bp->sliceloc == sloc)  // already have counted this location
	  break;
	if(sloc < bp->sliceloc) {
	  if((newbean = (struct bean *)malloc(sizeof(struct bean))) == NULL) {
	    perror("malloc failed");
	    exit(1);
	  }
	  if(bp->prev)
	    bp->prev->next = newbean;
	  newbean->prev = bp->prev;
	  newbean->next = bp;
	  bp->prev = newbean;
	  newbean->sliceloc = sloc;
	  if(bp == beanstar)
	    beanstar = newbean;
	  break;
	}
	if(!bp->next) {  // new node at end of chain
	  if((newbean = (struct bean *)malloc(sizeof(struct bean))) == NULL) {
	    perror("malloc failed");
	    exit(1);
	  }
	  newbean->next = (struct bean *)0;
	  newbean->prev = bp;
	  bp->next = newbean;
	  newbean->sliceloc = sloc;
	  break;
	}
      }
    }
  }
  for(i=0, bp = beanstar; bp; bp = bp->next,i++) {
    bp->order = i;
  }
  nbeans = i;
  for(ip1 = basep; ip1; ip1 = ip1->next)
    ip1->done = 0;
}


unlink(ip)
     struct image *ip;
{
  if(ip == basep)
    basep = ip->next;
  if(ip->prev)
    ip->prev->next = ip->next;
  if(ip->next)
    ip->next->prev = ip->prev;

  ip->next = (struct image *)0;
  ip->prev = (struct image *)0;
}


deleteim(ip)
     struct image *ip;

{
  if(ip == basep)
    basep = ip->next;

  if(ip->next)
    ip->next->prev = ip->prev;
  if(ip->prev)
    ip->prev->next = ip->next;

  if(ip->filename)
    free(ip->filename);
  if(ip->sliceloc)
    free(ip->sliceloc);
  if(ip->imagepos)
    free(ip->imagepos);
  if(ip->pixelspacing)
    free(ip->pixelspacing);
  if(ip->slicethickness)
    free(ip->slicethickness);
  if(ip->acqtime)
    free(ip->acqtime);
  if(ip->mediasopUID)
    free(ip->mediasopUID);
  if(ip->sopUID)
    free(ip->sopUID);
  if(ip->seriesUID)
    free(ip->seriesUID);
  if(ip->image)
    free(ip->image);
  if(ip->real)
    free(ip->real);
  if(ip->imag)
    free(ip->imag);
  if(ip->b0)
    free(ip->b0);
  if(ip->water)
    free(ip->water);
  if(ip->fat)
    free(ip->fat);
  if(ip->name)
    free(ip->name);
  free(ip);
}

doitman(index)
     int index;

{
  struct image *ip,*ip1,*ip2,*ip3,*ip4;
  struct image *ips[128];  // how many different echo times could we really have?
  int nips;
  char *thisloc;
  double tinc_wanted,diffmin;
  int i,j,k;
  double re,im,maxmag;
  double *ddumm,*ddump;
  double *mag,*phase;
  double denom,num,phaseoff;
  double ftsa,fatfac,waterfac;
  char fname[256];
  FILE *tfp;
  int nrows,ncols,resize;
  int firstrow,firstcol;
  double complex I1,I2,A,B,C,D,F,W;
  double *fatr,*fati,*waterr,*wateri;
  double round(),atof();
  struct bean *bp;
  int fwidth;
  double dphase,dtemp;
  double ddum;

  for(bp = beanstar; bp->next; bp = bp->next)
    if(bp->order == index)
      break;
  if(bp->order != index)
    fprintf(stderr,"could not find bean for requested index %d, bailing...\n",index);
  
  nips = 0;
  for(ip = basep;ip;ip = ip->next)
    if(!byimagenumber) {
      if(bp->sliceloc == atof(ip->sliceloc))
	break;
    }
    else {
      if(bp->sliceloc == ip->imageno)
	break;
    }
  if(!ip) {
    printf("could not find any images at slice location %g, bean index %d\n",bp->sliceloc,bp->order);
    return;
  }
  ips[nips++] = ip;
  for(ip1 = ip->next; ip1;ip1 = ip1->next) {
    if((ip1->rows != ip->rows) || (ip1->cols != ip->cols)
       || (ip1->accession != ip->accession) 
       || strcmp(ip1->pixelspacing,ip->pixelspacing)
       || strcmp(ip1->slicethickness,ip->slicethickness) || (ip1->tr != ip->tr)) {
      if(!strncmp(ip1->sliceloc,ip->sliceloc,4)) {
	printf("same locations mismatch at %s, filenames %s and %s\n",ip1->sliceloc,ip1->filename,ip->filename);
	printf("\tfile %s:\n\t\trows = %d\tcols = %d\n",ip1->filename,ip1->rows,ip1->cols);
	printf("\t\taccession = %d\tpixelspacing = %s\n",ip1->accession,ip1->pixelspacing);
	printf("\t\tslicethickness = %s\ttr = %g\n",ip1->slicethickness,ip1->tr);
	printf("\tfile %s:\n\t\trows = %d\tcols = %d\n",ip->filename,ip->rows,ip->cols);
	printf("\t\taccession = %d\tpixelspacing = %s\n",ip->accession,ip->pixelspacing);
	printf("\t\tslicethickness = %s\ttr = %g\n",ip->slicethickness,ip->tr);
	if(reallysloppy)
	  printf("(ignoring this mismatch)\n");
      }
      if(byimagenumber && (ip1->imageno == ip->imageno)) {
	printf("same image number mismatch at image number %d, filenames %s and %s\n",ip1->imageno,ip1->filename,ip->filename);
	printf("\tfile %s:\n\t\trows = %d\tcols = %d\n",ip1->filename,ip1->rows,ip1->cols);
	printf("\t\taccession = %d\tpixelspacing = %s\n",ip1->accession,ip1->pixelspacing);
	printf("\t\tslicethickness = %s\ttr = %g\n",ip1->slicethickness,ip1->tr);
	printf("\tfile %s:\n\t\trows = %d\tcols = %d\n",ip->filename,ip->rows,ip->cols);
	printf("\t\taccession = %d\tpixelspacing = %s\n",ip->accession,ip->pixelspacing);
	printf("\t\tslicethickness = %s\ttr = %g\n",ip->slicethickness,ip->tr);
	if(reallysloppy)
	  printf("(ignoring this mismatch)\n");
      }
      if(!reallysloppy)
	continue;
    }
    
    if(strcmp(ip1->sliceloc,ip->sliceloc) || strcmp(ip1->imagepos,ip->imagepos))
      if(!byimagenumber)
	continue;
      else
	{
	  if(ip1->imageno != ip->imageno)
	    continue;
	  else
	    printf("matching image number %d for slice locations %s and %s\n\t(imagepos = %s vs %s)\n",ip->imageno,ip->sliceloc,ip1->sliceloc,ip->imagepos,ip1->imagepos);
	}
    if(ip1->larmor_freq != ip->larmor_freq)
      if(fix) {
	printf("attempting phase correction of frequency mismatch of %g Hz\n",(ip1->larmor_freq - ip->larmor_freq)*1.0e6);
	dphase = (ip1->larmor_freq - ip->larmor_freq)*1.0e6 * ip1->te * 1.0e-3 * 2.0 * M_PI;
	for(i=0; i < ip1->rows * ip1->cols; i++) {
	  dtemp = ip1->real[i];
	  ip1->real[i] = ip1->real[i] * cos(dphase) + ip1->imag[i] * sin(dphase);
	  ip1->imag[i] = ip1->imag[i] * cos(dphase) - dtemp * sin(dphase);
	}
	ip1->larmor_freq = ip->larmor_freq;
      }
      else {
	if(!sloppy && !reallysloppy) {
	  printf("same locations frequency mismatch of %g Hz at %s, filenames %s and %s\n",(ip1->larmor_freq - ip->larmor_freq)*1.0e6,ip1->sliceloc,ip1->filename,ip->filename);
	  continue;
	}
	else
	  printf("(center freq mismatch %.10f vs %.10f (delta %g Hz))\n",ip1->larmor_freq,ip->larmor_freq,(ip1->larmor_freq - ip->larmor_freq)*1.0e6);
      }
    ips[nips++] = ip1;
  }
  for(i=0; i<nips; i++)
    unlink(ips[i]);
  if(!byimagenumber)
    printf("%d images for location %s\n",nips,ips[0]->sliceloc);
  else
    printf("%d images for image number %d\n",nips,ips[0]->imageno);

  if(nips == 1) {
    if(!byimagenumber)
      printf("only one image for location %s, skipping...\n",ips[0]->sliceloc);
    else
      printf("only one image for image number %d, skipping...\n",ips[0]->imageno);
    if(bp->next)
      bp->next->prev = bp->prev;
    if(bp->prev)
      bp->prev->next = bp->next;
    if(bp == beanstar)
      beanstar = bp->next;
    free(bp);
    return;
  }
  tinc_wanted = 1000.0 / (ips[0]->larmor_freq * fabs(fwdiff));  // larmor_freq is in MHz
  diffmin = tinc_wanted;
  
  ip1 = ip2 = (struct image *)0;  // find appropriate series for b0 map
  for(i = 0; i < nips; i++)
    for(j = i+1; j < nips; j++) {
      if(ips[i]->te == ips[j]->te) {
	printf("averaging series %d image %d with series %d image %d\n",ips[i]->series,ips[i]->imageno,ips[j]->series,ips[j]->imageno);
	for(k=0; k < ips[i]->rows * ips[i]->cols; k++) {
	  ips[i]->real[k] = (ips[i]->real[k] + ips[j]->real[k]) / 2.0;
	  ips[i]->imag[k] = (ips[i]->imag[k] + ips[j]->imag[k]) / 2.0;
	}
	for(k=j; k<(nips-1); k++)
	  ips[k] = ips[k+1];
	nips--;
	j--;
	continue;
      }
      if(fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted) < (diffmin + SLOP)) {
	if(!ip1) { // no competition yet
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip1 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip2 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	  continue;
	}
	if((fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted) + SLOP) < diffmin) { // clear winner
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip1 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip2 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	  continue;
	}
	// would rather use in-phase echo time if possible
	if(min(fabs(fmod(ips[i]->te,tinc_wanted)),fabs(tinc_wanted - fmod(ips[i]->te,tinc_wanted)))
	   < min(fabs(fmod(ip1->te,tinc_wanted)),fabs(tinc_wanted - fmod(ip1->te,tinc_wanted)))) {
	  // this echo time is closer to in-phase and about as good time difference, use it
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip1 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip2 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	}
      }
    }
  if(!byimagenumber)
    printf("location %s, echo times %g and %g, tinc = %g, ideal tinc = %g\n",ip1->sliceloc,ip1->te,ip2->te,ip2->te - ip1->te,tinc_wanted);
  else
    printf("image number %d, echo times %g and %g, tinc = %g, ideal tinc = %g\n",ip1->imageno,ip1->te,ip2->te,ip2->te - ip1->te,tinc_wanted);
  maxmag = 0.0;
  
  if((mag = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((phase = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }

  for(i=0; i<ip1->rows*ip1->cols; i++) {
    //    mag[i] = hypot(ip1->real[i],ip1->imag[i]) * hypot(ip2->real[i],ip2->imag[i]);
    mag[i] = sqrt(ip1->real[i]*ip1->real[i]+ip1->imag[i]*ip1->imag[i]+ip2->real[i]*ip2->real[i]+ip2->imag[i]*ip2->imag[i]);
    if(mag[i] > maxmag)
      maxmag = mag[i];
    // not normalized by amplitude because angle is all we care about
    re = ip1->real[i] * ip2->real[i] + ip1->imag[i] * ip2->imag[i];
    im = ip1->real[i] * ip2->imag[i] - ip1->imag[i] * ip2->real[i];
    phase[i] = atan2(im,re);
  }
  for(i=0; i<ip1->rows * ip1->cols; i++) {  // make sure shorts in punwrap don't overflow.
    mag[i] /= maxmag;                 // believe me, you don't want your shorts to overflow.
    mag[i] *= 1000.0;
  }
  
  printf("phase unwrapping...");
  fflush(stdout);
  
  nrows = ncols = 0;   // resize to power of two for phase unwrapping
  resize = 0;
  if(!ispow2(ip1->rows)) {
    nrows = (int)ceil(log((double)ip1->rows) / log(2.0));
    nrows = (int)pow(2.0,(double)nrows);
  }
  if(!ispow2(ip1->cols)) {
    ncols = (int)ceil(log((double)ip1->cols) / log(2.0));
    ncols = (int)pow(2.0,(double)ncols);
  }
  if(ncols || nrows) {
    resize = 1;
    if(!ncols)
      ncols = ip1->cols;
    if(!nrows)
      nrows = ip1->rows;
    if((ddumm = (double *)malloc(sizeof(double)*nrows*ncols)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    if((ddump = (double *)malloc(sizeof(double)*nrows*ncols)) == NULL) {
      perror("malloc failed");
      exit(1);
    }
    // could do Fourier resizing with routine from libglen.a
    //	resize(mag,phase,ip1->rows,ip1->cols,ddumm,ddump,nrows,ncols);
    // instead, just pad with zeros, much faster
    
    
    firstrow = (nrows - ip1->rows)/2 + (nrows - ip1->rows)%2;
    firstcol = (ncols - ip1->cols)/2 + (ncols - ip1->cols)%2;
    
    for(j=0; j<nrows; j++) {
      for(k=0; k<ncols; k++) {
	if((j >= firstrow) && (j < firstrow + ip1->rows) && (k >= firstcol) && (k < firstcol + ip1->cols)){
	  ddumm[j*ncols+k] = mag[(j - firstrow)*ip1->cols + (k - firstcol)];
	  ddump[j*ncols+k] = phase[(j - firstrow)*ip1->cols + (k - firstcol)];
	}
	else {
	  ddumm[j*ncols+k] = 0.0;
	  ddump[j*ncols+k] = 0.0;
	}
      }
    }
  }
  else {  // no resizing needed
    ddump = phase;
    ddumm = mag;
    ncols = ip1->cols;
    nrows = ip1->rows;
  }

  punwrapstub2(ddump,ddumm,nrows,ncols,1,10.0 * threshlow,10.0 * threshhigh);
  
  if(resize) {
    for(j=0; j<ip1->rows; j++)
      for(k=0; k<ip1->cols; k++)
	phase[j*ip1->cols+k] = ddump[(firstrow + j)*ncols + (firstcol + k)];
    free(ddump);
    free(ddumm);
  }
  for(i=0; i<ip1->rows*ip1->cols; i++) {
    if(threshmask == -1.0)
      mag[i] = (mag[i] > (10.0 * threshlow))?(1000.0):(0.0);
    else
      mag[i] = (mag[i] > (10.0 * threshmask))?(1000.0):(0.0);
  }
  // mag is now mask
  
  //  sprintf(fname,"Mask_");
  //  writeDicom(fname,ip1,mag,firstseries+5,ip1->imageno,"Mask",poff,cmdline);

  printf(" done\n");
  
  // phase is currently in units of radians.  Convert to parts per billion, i.e. 1000 = 1ppm.
  num = denom = 0.0;
  for(i=0; i<ip1->rows*ip1->cols; i++)
    if(mag[i]) {
      num += phase[i];
      denom += 1.0;
    }
    else
      phase[i] = 0.0;
  
  // average phase in non-noise pixels of b0 map.
  
  if(!denom)
    bp->b0mean = 0.0;
  else
    bp->b0mean = num / denom;
  //  printf("slice location %g, average of b0 map = %g radians.\n",bp->sliceloc,bp->b0mean);
  if(firstb0done) {  // compare to adjacent slice to limit jumps in b0
    if(index > nbeans/2) { // look back
      //      printf("\tcomparing to previous slice of location %g, %g radians.\n",bp->prev->sliceloc,bp->prev->b0mean);
      phaseoff = -2.0 * M_PI * round((bp->b0mean - bp->prev->b0mean) / (2.0 * M_PI));
      //      printf("\tadding phase angle %g to slice location %g\n",phaseoff,bp->sliceloc);
    }
    else { // look forward
      //      printf("\tcomparing to next slice of location %g, %g radians.\n",bp->next->sliceloc,bp->next->b0mean);
      phaseoff = -2.0 * M_PI * round((bp->b0mean - bp->next->b0mean) / (2.0 * M_PI));
      //      printf("\adding phase angle %g to slice location %g\n",phaseoff,bp->sliceloc);
    }
  }
  else { // potential first slice, make in range -pi to pi
    if(fabs(bp->b0mean) > M_PI)
      phaseoff = -2.0 * M_PI * round(bp->b0mean / (2.0 * M_PI));
    if(oops)
      phaseoff += 2.0 * M_PI * (double)oops;
  }
  if(phaseoff)
    for(i=0; i<ip1->rows*ip1->cols; i++)
      if(mag[i])
	phase[i] += phaseoff;
  
  // check, remove this later, replace with bp->b0mean += phaseoff
  
  num = denom = 0.0;
  for(i=0; i<ip1->rows*ip1->cols; i++)
    if(mag[i]) {
      num += phase[i];
      denom += 1.0;
    }
  if(!denom)
    bp->b0mean = 0.0;
  else
    bp->b0mean = num / denom;
  //  printf("slice location %g, now has average of b0 map = %g radians.\n",bp->sliceloc,bp->b0mean);
  
  // convert to units of parts per billion
  
  for(i=0; i<ip1->rows*ip1->cols; i++)
    phase[i] = phase[i] / ((ip2->te - ip1->te) * 1.0e-3) / (2.0 * M_PI) / ip1->larmor_freq * 1000.0;
  
  ip1->b0 = phase;
  //  fwrite(ip1->b0,sizeof(double),ip1->rows*ip1->cols,debugfp);
  if(dob0) {
    if(rawimages) {
      sprintf(fname,"rawB0\0");
      writeRaw(fname,ip1,ip1->b0);
    }
    
    if(dicomimages) {
      sprintf(fname,"B0\0");
      writeDicom(fname,ip1,ip1->b0,firstseries,ip1->imageno,"B0 map",poff,cmdline); // series 100 start at image 1
    }
  }
  // Do three-point fat-water decomposition
  
  if(nips < 3) {
    if(!byimagenumber)
      printf("not enough echo times to do Dixon recon for location %.5s, skipping...\n",ip1->sliceloc);
    else
      printf("not enough echo times to do Dixon recon for image number %d, skipping...\n",ip1->imageno);
    if(bp->next)
      bp->next->prev = bp->prev;
    if(bp->prev)
      bp->prev->next = bp->next;
    if(bp == beanstar)
      beanstar = bp->next;
    free(bp);
    return;
  }
  
  // pick two likely echoes
  
  tinc_wanted = 1000.0 / (ips[0]->larmor_freq * fabs(fwdiff)) / 2.0;  // larmor_freq is in MHz
  diffmin = tinc_wanted;
  ip3 = ip4 = (struct image *)0;
  for(i = 0; i < nips; i++)
    for(j = i+1; j < nips; j++) {
      if(ips[i] == ip1)
	if(ips[j] == ip2)  // can't use these two again, not independent.
	  continue;
      if(ips[i] == ip2)    // ditto
	if(ips[j] == ip1)
	  continue;
      
      //      printf("ips[%d]->te = %g, ips[%d]->te = %g, diffmin = %g\n",i,ips[i]->te,j,ips[j]->te,diffmin);
      if(fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted) < (diffmin + SLOP)) {
	if(!ip3) {
	  // no competition
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip3 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip4 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	  //	  printf("choosing echo times %g and %g as first comers\n",ip3->te,ip4->te);
	  continue;
	}
	if((fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted) + SLOP) < diffmin) { // clear winner
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip3 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip4 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	  //	  printf("choosing echo times %g and %g as better\n",ip3->te,ip4->te);
	  continue;
	}
	// take earliest echos possible
	if(min(ips[i]->te,ips[j]->te) < min(ip3->te,ip4->te)) {
	  diffmin = fabs(fabs(ips[i]->te - ips[j]->te) - tinc_wanted);
	  ip3 = (ips[i]->te > ips[j]->te)?(ips[j]):(ips[i]);
	  ip4 = (ips[i]->te > ips[j]->te)?(ips[i]):(ips[j]);
	  //	  printf("choosing echo times %g and %g as just as good and earlier\n",ip3->te,ip4->te);
	}	    
      }
    }
  if(!pixelx) {
    printf("pixel spacing: %s\n",ip3->pixelspacing);
    if(sscanf(ip3->pixelspacing,"%lf\\%lf",&pixelx,&pixely) != 2) {
      perror("sscanf error");
      exit(1);
    }
  }
  printf("fw decomp te = %g and te = %g, tinc = %g, ideal = %g\n",ip3->te,ip4->te,ip4->te - ip3->te,tinc_wanted);
  if(poff)
    printf("adding %g degrees of phase to image with echo time %g\n",poff,ip4->te);
  if((fatr = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((waterr = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((fati = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  if((wateri = (double *)malloc(sizeof(double)*ip1->cols*ip1->rows)) == NULL) {
    perror("malloc failed");
    exit(1);
  }
  for(i = 0; i < ip3->rows * ip3->cols; i++) {
    // I never knew you could do complex math with newer c99 compilers.  Cool.
    I1 = ip3->real[i] + I*ip3->imag[i];
    I2 = ip4->real[i] + I*ip4->imag[i];

    if(poff)
      I2 = creal(I2)*cos(poff*M_PI/180) - cimag(I2)*sin(poff*M_PI/180) + 
	I*(creal(I2)*sin(poff*M_PI/180) + cimag(I2)*cos(poff*M_PI/180));

    A = cexp(I * 2.0 * M_PI * fwdiff * ip3->te * ip3->larmor_freq * 1.0e-3);
    B = cexp(I * 2.0 * M_PI * fwdiff * ip4->te * ip3->larmor_freq * 1.0e-3);
    C = cexp(I * -2.0 * M_PI * ip1->b0[i]/1000.0 * ip3->te * ip3->larmor_freq * 1.0e-3);
    D = cexp(I * -2.0 * M_PI * ip1->b0[i]/1000.0 * ip4->te * ip3->larmor_freq * 1.0e-3);
    W = (I1 * A * C - I2 * B * D) / (A - B);
    F = I1*A*C - W*A;

    //    convert to mag and phase
    waterr[i] = cabs(W);
    // follow Siemens phase image convention of 0-4096 <--> 0-2pi
    wateri[i] = fmod(carg(W) + 2.0 * M_PI, 2.0 * M_PI) * 4096.0 / (2.0 * M_PI);
    fatr[i] = cabs(F);
    fati[i] = fmod(carg(F) + 2.0 * M_PI, 2.0 * M_PI) * 4096.0 / (2.0 * M_PI);
  }

  /*  fwrite(waterr,sizeof(double),ip3->rows*ip3->cols,debugfp);
  fwrite(wateri,sizeof(double),ip3->rows*ip3->cols,debugfp);
  fwrite(fatr,sizeof(double),ip3->rows*ip3->cols,debugfp);
  fwrite(fati,sizeof(double),ip3->rows*ip3->cols,debugfp);
  printf("rows = %d, cols = %d, slice loc = %g\n",ip3->rows,ip3->cols,bp->sliceloc);
  */
  if(dowater) {
    if(rawimages) {
      sprintf(fname,"rawWater_%03d\0",ip3->imageno);
      //      fwrite(waterr,sizeof(double),ip3->rows*ip3->cols,debugfp);
      writeRaw(fname,ip3,waterr);
    }
    if(dicomimages) {
      sprintf(fname,"Water_\0");
      writeDicom(fname,ip3,waterr,firstseries+1,ip3->imageno,"Water mag",poff,cmdline);  // series 101
    }
    if(dophase) {
      if(rawimages) {
	sprintf(fname,"rawWater_phase_%03d\0",ip3->imageno);
	writeRaw(fname,ip3,wateri);
      }
      if(dicomimages){
	sprintf(fname,"Water_phase_\0");
	writeDicom(fname,ip3,wateri,firstseries+2,ip3->imageno,"Water phase",poff,cmdline);
      }
    }
  }
  if(dofat) {
    if(rawimages) {
      sprintf(fname,"rawFat_%03d\0",ip3->imageno);
      writeRaw(fname,ip3,fatr);
    }
    if(dicomimages){
      sprintf(fname,"Fat_\0");
      writeDicom(fname,ip3,fatr,firstseries+3,ip3->imageno,"Fat mag",poff,cmdline);  // series 103
    }
    if(dophase) {
      if(rawimages) {
	sprintf(fname,"rawFat_phase_%03d\0",ip3->imageno);
	writeRaw(fname,ip3,fati);
      }
      if(dicomimages) {
	sprintf(fname,"Fat_phase_\0");
	writeDicom(fname,ip3,fati,firstseries+4,ip3->imageno,"Fat phase",poff,cmdline);
      }
    }
  }
  // raw format for Scion Image
  /*  fwidth = min(strlen(ip3->sliceloc),5);
  printf("strlen = %d, str = %s...\n",strlen(ip3->sliceloc),ip3->sliceloc);
  sprintf(fname,"Fat_raw_%.*s\0",fwidth,ip3->sliceloc);
  writeRaw(fname,ip3,fatr);
  */

  // Calculation of fat-water fraction
  fatfac = (1.0 - cos(ip3->flipangle*M_PI/180.0)*exp(-ip3->tr/t1fat))
    / (sin(ip3->flipangle*M_PI/180.0)*(1.0 - exp(-ip3->tr/t1fat)));
  waterfac = (1.0 - cos(ip3->flipangle*M_PI/180.0)*exp(-ip3->tr/t1wat))
    / (sin(ip3->flipangle*M_PI/180.0)*(1.0 - exp(-ip3->tr/t1wat)));

  for(i=0; i<ip3->rows*ip3->cols; i++) {
    fatr[i] = fatr[i]*fatfac;
    waterr[i] = waterr[i]*waterfac;
    ftsa = fatr[i] / (fatr[i] + waterr[i]);
    fatr[i] = 1000.0 * ftsa / (1.138 - 0.339*ftsa);  
    // using fati for water volume fraction
    fati[i] = 1000.0 - fatr[i];
    if(fati[i] < 0.0)
      fati[i] = 0.0;
    if(!mag[i]) {
      fatr[i] = 0.0;
      fati[i] = 0.0;
    }
    // fat volume fraction in units of 1/10 of a percent
  }
  bp->nvalidpts = 0;
  bp->meanfwfrac = 0.0;
  bp->rows = ip3->rows;
  bp->cols = ip3->cols;
  for(i=0; i<ip3->rows*ip3->cols; i++) {
    if(mag[i]) {
      bp->nvalidpts++;
      bp->meanfwfrac += fatr[i];
    }
  }
  if(bp->nvalidpts != 0)
    bp->meanfwfrac /= (double)bp->nvalidpts;

  if(dofwfrac) {
    if(rawimages) {
      sprintf(fname,"rawFat_frac_%03d\0",ip3->imageno);
      writeRaw(fname,ip3,fatr);
    }
    if(dicomimages) {
      sprintf(fname,"Fat_frac_\0");
      writeDicom(fname,ip3,fatr,firstseries+5,ip3->imageno,"Fat fraction",poff,cmdline);  // series 105
    }
  }

  if(dowatfrac) {
    if(rawimages) {
      sprintf(fname,"rawWat_frac_%03d\0",ip3->imageno);
      writeRaw(fname,ip3,fati);
    }
    if(dicomimages) {
      sprintf(fname,"Wat_frac_\0");
      writeDicom(fname,ip3,fati,firstseries+7,ip3->imageno,"Water fraction",poff,cmdline);  // series 107
    }
  }

  if(!firstb0done)
    firstb0done++;
  
  for(i=0; i<nips; i++)
    deleteim(ips[i]);
  
  free(mag);
  free(waterr);
  free(wateri);
  free(fatr);
  free(fati);
}


indent()
{
  int i;

  for(i=0; i<level; i++)
    printf("\t");
}
