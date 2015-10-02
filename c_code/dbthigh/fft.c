/*
 * file name	: fft.c
 *
 * Author	: S. M. Song
 *
 * Description	: This file contains subroutines for 1-D variable
 *			 size FFT via decimation-in-time in-place algorithm.
 *
 * Copyright (C) 1993, Samuel M. Song.  All rights reserved.  
 *
 */
#include<stdlib.h>

#include	<math.h>
#include	<stdio.h>
#include	<stddef.h>
#include	"fft.h"

/*
	Below externals are initialized by a call to fftinit() automatically.
	Remember to call fftexit() to clean up.
*/

COMPLEX		*Ffftc,*Fintpe,*Fintpo;
short		*Fbitrev;
int		sizeo2,lg2sze;
float		gain;

/*
 * Subroutine	: fft1d
 *
 * Description	: This subroutine performs 1-D variable size FFT via
 *			decimation-in-time in-place algorithm.
 *		  The output overwrites the input array.
 *
 */
int fft1d(ainp)
	COMPLEX		*ainp;

{
	short		ipass;
	COMPLEX		prod;
	COMPLEX		*pintin,*pintout,*pinttemp,*pout;
	register short	cmask, i, j;

#define		CMULT(A,B,C)	A.real = B.real*C.real - B.imag*C.imag, \
				A.imag = B.real*C.imag + B.imag*C.real
#define		CADD(A,B,C)	A.real = B.real + C.real, \
				A.imag = B.imag + C.imag
#define		CSUB(A,B,C)	A.real = B.real - C.real, \
				A.imag = B.imag - C.imag

 						/* FFT pass 1 */
	for (i = 0, j = sizeo2; i < sizeo2; i++, j++) {
		CADD(Fintpo[2*i],ainp[i],ainp[j]);
		CSUB(Fintpo[2*i+1],ainp[i],ainp[j]);
	}

	pintin = Fintpo; pintout = Fintpe;
	cmask = 0;				/* the rest */
	for (ipass = 1; ipass < lg2sze; ipass++) {
		cmask += (1 << (ipass-1));
		for (i = 0, j = sizeo2; i < sizeo2; i++, j++) {
			CMULT(prod,pintin[j],Ffftc[i&cmask]);
			CADD(pintout[2*i],pintin[i],prod);
			CSUB(pintout[2*i+1],pintin[i],prod);
		}
		pinttemp = pintout;		/* rotate (swap) buffers */
		pintout  = pintin;
		pintin = pinttemp;
	}

/*
	Bit reverse the data and put 1/sqrt(size) on the output.
	If number of passes were even, then output is in array Fintpe;
	else (odd number of passes), the output is in array Fintpo  
*/

	if ((lg2sze%2) == 0) pout = Fintpe;		/* select buffer */
	else pout = Fintpo;

	for (i = 0; i < sizeo2; i++) {
		ainp[2*Fbitrev[i]].real = gain * pout[i].real;
		ainp[2*Fbitrev[i]].imag = gain * pout[i].imag;
		ainp[2*Fbitrev[i]+1].real = gain * pout[i+sizeo2].real;
		ainp[2*Fbitrev[i]+1].imag = gain * pout[i+sizeo2].imag;
	}
}

/*
 * Subroutine	: fftinit
 *
 * Description	: This subroutine initializes various arrays and variables
 *			for fft1d routine.
 *
 */

int fftinit(xdir,size)

int	xdir,size;			/* SMS intflg was real */

{
	register int	i,j;
	short		length = 1;
	float		pi = 4. * atan(1.),twopi = 2. * pi,phase;

        sizeo2 = size / 2;			/* external variables */
        lg2sze = (int) (log((double) (size)) / log(2.0));

	if (xdir == 1)
		gain = 1.;
	else gain = 1. / (float) (size);

/*
	allocate arrays
*/

	Fbitrev = (short *) (calloc((size_t) (size),
					(size_t) (sizeof(*Fbitrev))));
	Ffftc  = (COMPLEX *) (calloc((size_t) (sizeo2),
					(size_t) (sizeof(*Ffftc))));
	Fintpe = (COMPLEX *) (calloc((size_t) (size),
					(size_t) (sizeof(*Fintpe))));
	Fintpo = (COMPLEX *) (calloc((size_t) (size),
					(size_t) (sizeof(*Fintpo))));

	if ((Fbitrev && Ffftc && Fintpe && Fintpo) == 0) {
		(void) fprintf(stderr,"fftinit(): cannot allocate arrays\n");
		(void) exit(-1);
	}

/*
	Build a table of bit reversed numbers
*/

        Fbitrev[0] = 0;
	for (i = 0; i < lg2sze-1; i++) {
		for (j = 0; j < length; j++) {
			Fbitrev[j] = 2*Fbitrev[j];
			Fbitrev[j+length] = Fbitrev[j] + 1;
		}
		length = 2*length;
	}

/*
	Build the table of FFT coeff (complex phasors) in bit rev format
*/

	for (i = 0; i < sizeo2; i++) {
		phase = twopi * ((float) (Fbitrev[i]))
				/ ((float) (size));
		phase = -xdir * phase;
		Ffftc[i].real = cos((double) (phase));
		Ffftc[i].imag = sin((double) (phase));
	}
}

/*
 * Subroutine	: fftexit
 *
 * Description	: This subroutine cleans-up the mess
 *
 */

int	fftexit()
{
	free(Fbitrev);
	free(Ffftc);
	free(Fintpo);
	free(Fintpe);
}
