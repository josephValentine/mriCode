/*
 * file name	: cost.c
 *
 * Author	: S. M. Song
 *
 * Description	: This file contains subroutines for variable size
 *			cosine transform algorithm via the FFT.
 * 
 * Copyright (C) 1993, Samuel M. Song.  All rights reserved. 
 *
 */
#include<stdlib.h>
#include	<math.h>
#include	<stdio.h>
#include	<fcntl.h>
#include	"fft.h"

/*
	Below external is initialized by a call to costinit(size) automatically.
	Remember to call costexit() to clean up.
*/

COMPLEX		*Cphasef;
int		Cold_size = 0;

/*
 * Subroutine	: cost3d
 *
 * Description	: This subroutine performs 3-D variable size cosine
 *			transform via the FFT.
 *
 */
int cost(ainp,sizex,sizey,sizez,xdir)
	float		*ainp;
	int		sizex,sizey,sizez,xdir;
{
	register int	i,j,k;
	COMPLEX		*temp;

	for (k = 0; k < sizez; k++) {			/* in x and y-dir */
		cost2d(ainp+k*(sizex*sizey),sizex,sizey,xdir);
	}

	if (sizez != 1) {
		if ((temp = (COMPLEX *) (calloc((size_t) (2*sizez),
					(size_t) (sizeof(*temp))))) == NULL) {
			(void) fprintf(stderr,
				"cost3d(): cannot allocate arrays\n");
			(void) exit(-1);
		}
		for (i = 0; i < sizex; i++) {		/* in z-direction */
		for (j = 0; j < sizey; j++) {
			for (k = 0; k < sizez; k++) {
				temp[k].real =
					 ainp[k*sizex*sizey + j*sizex + i];
				temp[k].imag = 0.;
				temp[k+sizez].real = 0.;	
				temp[k+sizez].imag = 0.;	
			}
			cost1d(temp,sizez,xdir);
			for (k = 0; k < sizez; k++) {
				ainp[k*sizex*sizey + j*sizex + i] = 
					temp[k].real;
			}
		}
		}
		free(temp);
	}

}
/*
 * Subroutine	: cost2d
 *
 * Description	: This subroutine performs 2-D variable size cosine
 *			transform via the FFT.
 *
 */
int cost2d(ainp,sizex,sizey,xdir)
	float		*ainp;
	int		sizex,sizey,xdir;
{
	register int	i,j;
	int		offset,maxsize;
	COMPLEX		*temp;
	//	FILE *debug;

	maxsize = (sizex > sizey) ? sizex : sizey;
	temp = (COMPLEX *) (calloc((size_t) (2 * maxsize),
					(size_t) (sizeof(*temp))));

	if (temp == NULL) {
		(void) fprintf(stderr,"cost2d(): cannot allocate array\n");
		(void) exit(-1);
	}

	for  (j = 0; j < sizey; j++) {		/* xform on rows */
		offset = sizex * j;
		for (i = 0; i < sizex; i++) {
			temp[i].real = ainp[i+offset];
			temp[i].imag = 0.;
			temp[i+sizex].real = 0.;
			temp[i+sizex].imag = 0.;
		}
		cost1d(temp,sizex,xdir);
		for (i = 0; i < sizex; i++) {
			ainp[i+offset] = temp[i].real;
		}
	}

	for  (i = 0; i < sizex; i++) {		/* on columns */
		for (j = 0; j < sizey; j++) {
			offset = sizex * j;
			temp[j].real = ainp[i+offset];
			temp[j].imag = 0.;
			temp[j+sizey].real = 0.;
			temp[j+sizey].imag = 0.;
		}
		cost1d(temp,sizey,xdir);
		for (j = 0; j < sizey; j++) {
			offset = sizex * j;
			ainp[i+offset] = temp[j].real;
		}
	}

	free(temp);
	//	debug = fopen("./debug","w");
	//	fwrite(ainp,sizeof(float),sizex*sizey*2,debug);
	//	exit(0);
}
/*
 * Subroutine	: cost1d
 *
 * Description	: This subroutine performs 1-D variable size cosine
 *			transform via the FFT.
 *		  The output overwrites the input array.
 *
 *		  Note well:
 *			input and output arrays are COMPLEX twice the length
 *
 */
int cost1d(cinp,size,xdir)
	COMPLEX		*cinp;
	int		size,xdir;
{
	float		sqrt2 = sqrt(2.);
	float		sqrt2on = sqrt(2./(float)(size));
	register int	i;

	if (Cold_size != size) {	/* if size changed	  */
		if (Cold_size != 0)	/*	if not first time */
			costexit();	/*		clean up  */
		costinit(size);		/*	initialize	*/
		Cold_size = size;	/*	reset size	*/
	}
	
	if (xdir == 1) {
	  fft1d(cinp);
	  for (i = 0; i < size; i++) {
	    cinp[i].real = sqrt2on *
	      ( Cphasef[i].real * cinp[i].real
		- Cphasef[i].imag * cinp[i].imag);
	  }		/* sqrt(2/n) * real(phasef .* X(1:n)); */
	  cinp[0].real /= sqrt2;
	}
	else {
	  cinp[0].real /= sqrt2;
	  for (i = 0; i < size; i++) {
	    cinp[i].imag =  cinp[i].real * Cphasef[i].imag;
	    cinp[i].real = 	cinp[i].real * Cphasef[i].real;
	  }		/* X(1:n) .* phasef;	loop */
	  fft1d(cinp);
	  for (i = 0; i < size; i++) {
	    cinp[i].real = sqrt2on * cinp[i].real;
	  }		/* X = sqrt(2/n) * X(1:n); */
	}
}
/*
 * Subroutine	: costinit
 *
 * Description	: This subroutine initializes a phase factor array
 *			for the cosine transform.
 *
 */
int costinit(size)
	int		size;
{
	register int	i;
	float		mpi = -4. * atan(1.);
	float		phase;

	fftinit(1,2*size);		/* init FFT mods (forward FFTs only) */

	if ((Cphasef = (COMPLEX *) (calloc((size_t) (size),
				(size_t) (sizeof(*Cphasef))))) == NULL) {
		fprintf(stderr, "costinit(): cannot allocate arrays\n");
		exit(-1);
	}
	for (i = 0; i < size; i++) {
		phase = mpi * (float) (i) / ((float) (2*size));
		Cphasef[i].real = cos(phase);
		Cphasef[i].imag = sin(phase);
	
	}
}
/*
 * Subroutine	: costexit
 *
 * Description	: This subroutine cleans up the mess.
 *
 */
int costexit()
{
	Cold_size = 0;
	fftexit();
	free(Cphasef);
}
