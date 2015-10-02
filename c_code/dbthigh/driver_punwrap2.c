#include<stdlib.h>
#include 	<stdio.h>
#include 	<stddef.h>
#include	<fcntl.h>
#include	<math.h>
#include	<sys/types.h>
#include	<sys/stat.h>
#include	"punwrap.h"

#define		MAX_FILE_NAME	128
#define		DEBUG_IO	0

float	Rwval = 3141.0;			/* -1000 pi to 1000 pi */
int	Runweight = 0;
int	Rniter = 32;			/* number of iterations */
float	Rgnull = 0.0;			/* gradient nulling value */
int	Rthr_type = 0;			/* threshold type */
float	Rthrv1 = 3.0, Rthrv2 = 30.0;	/* threshold start and end */
int	R3d = 0;			/* process 3-D at once */

char	*Rfp1 = "I", *Rfp2 = "I";	/* input phase & mag file prefix */
char	*Rfm1 = "I", *Rfm2 = "I";
char	*Rfout = "J";			/* output file prefix */
int	Rimgtype;			/* file type 0(flat), 4(4x), 5(5x) */
long	Rfile_offset;			/* header offset */
int	Rfnametype;			/* file name type 0(I.001) 1(I.1) */

int	Rsp1 = 0, Rsp2 = 0,		/* file ext for phases and mag */
	Rsm1 = 0, Rsm2 = 0;

int	Rskip = 1;			/* file ext skips for phase and mag */

long	Rmonitor = 0;			/* monitoring flag: see punwrap.c */

IMAGE	*Ry;				/* unwrapped phase */
					/* program name, argv[0], sets this */
int	Rcoll_type = 0;

punwrapstub2(phaseim,magim,ydim,xdim,threshtype,threshl,threshh)
double *phaseim,*magim;
int ydim,xdim;
int threshtype;
float threshl,threshh;

{
	IMAGE	*x, *mag, *hdr;
	long	imgsize;		/* image size */
	long	lnxny;
	int	loop;			/* loop index */
	int	nimgs_pl;	/* loop count, n imgs read per loop */
	int	i;			/* misc */
	char	*readimg();
	char	*writeimg();
	int	curimgp1;
	float	wval_eff, scalef;
	float	gnull_eff;
	//	FILE *debug2;

	/* decode options */

	imgsize = xdim * ydim;
	lnxny   = (long) (xdim*ydim);

	x = (IMAGE *) calloc((size_t) (imgsize), (size_t) (sizeof(*x)));
	Ry = (IMAGE *) calloc((size_t) (imgsize), (size_t) (sizeof(*Ry)));
	mag = (IMAGE *) calloc((size_t) (imgsize), (size_t) (sizeof(*mag)));
	hdr = (IMAGE *) calloc((size_t) (Rfile_offset),
		(size_t) (sizeof(*hdr)));

	if (!(x && Ry && mag)) {
		fprintf(stderr, "driver_punwrap: cannot allocate all arrays\n");
		exit(-1);
	}

	for(i=0; i<imgsize; i++) {
	  mag[i] = (short)rint(magim[i]);
	  x[i] = (short)rint(phaseim[i] * 1000.0);
	  if(x[i] == -32768) 
	    printf("hey! x[%d] = %d, phaseim[%d] = %g\n",i,x[i],i,phaseim[i]);
	}
	//	debug2 = fopen("debug2","w");
	//	fwrite(x,sizeof(short),imgsize,debug2);
	//	fclose(debug2);
	wval_eff = Rwval;
	gnull_eff = Rgnull;
	
	//	punwrap(Ry, x, mag, wval_eff, xdim, ydim,
	//			1, Rniter, gnull_eff,
	//			Rthr_type, Rthrv1, Rthrv2, Rmonitor);
		punwrap(Ry, x, mag, wval_eff, xdim, ydim,
				1, Rniter, gnull_eff,
				threshtype, threshl, threshh, Rmonitor);
	for(i=0; i<imgsize; i++)
	  phaseim[i] = (double)Ry[i] / 1000.0;

	punwrapexit();				/* clean up */
	free(x), free(Ry), free(mag), free(hdr);

}


