#include<stdlib.h>
#include 	<stdio.h>
#include 	<stddef.h>
#include	<math.h>
#include	"punwrap.h"

/*
	Below externals initialized by a call to punwrapinit() automatically.
	User must call punwrapexit() to exit gracefully.
	See below.

*/

int     Pnx, Pny, Pnz, Pndim;		/* image dimensions */
long	Pnxny, Pnxnynz;
float	Pwval, Pwvalx2;
float	*Pthr;				/* ptr to thresholds */
FLOAT	*Pyls, *Pgrad, *Pcosx, *Pcosy, *Pcosz;

int	Pcntr = 0;			/* counter */

					/* monitor flags */
int	Catch_timing, Catch_thrs;
int	Catch_grad, Catch_div, Catch_neubc, Catch_yls;
int	Catch_magph, Catch_wraps;

punwrap(y, x, mag, wval, nx, ny, nz, niter, gnull, thr_type, thrv1, thrv2,
	monitor_flags)
	IMAGE	*y, *x, *mag;
	float	wval;
	int	nx, ny, nz, niter;
	float	gnull;
	int	thr_type;			/* zero OK */
	float	thrv1, thrv2;			/* zeros OK */
     long	monitor_flags;			/* zero OK */
{
  int	i;
  char	*writeimg();
  
  punwrapinit(wval, nx, ny, nz, niter, thr_type, thrv1, thrv2, mag);

  Pcntr++;
  
  for (i = 0; i < niter; i++)	/* update yls iteratively */
    lsphase(x,mag,Pyls,gnull,Pthr[i], i);

  /* adjust so that x = W[x] */
  padjust(y, Pyls, x, mag,		/* SMS threshold */
	  Pthr[0]);
}


/*
	Following macro used in lsphase() and padjust()
	
*/

#define	LOOP3D(Ind, Is,Ie,Js,Je,Ks,Ke,STMTS) {				\
		int Ind,J,K;						\
		for (K = Ks; K < Ke; K = K + Pnxny) {			\
			for (J = K + Js; J < K + Je; J = J + Pnx) {	\
				for (Ind = J + Is; Ind < J + Ie; Ind++) { \
			        	STMTS	\
				}	\
			}	\
		}	\
	}

/* 
 * Subroutine	: lsphase
 *
 * Description	: computes the least squares phase
 *
 */
lsphase(x, mag, yls, gnull, thr, ind)
	IMAGE	*x, *mag;
	FLOAT	*yls;
	float	gnull, thr;
	int	ind;
{
	int	i,is,ie,js,je,ks,ke;

#define	PWRAP(X)	((X) - Pwvalx2*floor(((X) + Pwval)/(Pwvalx2)))

/*
	compute the gradient: 	mag large -- wrap(grad(x))
				else	  -- grad(yls)
	This is a backward difference.  Compute all that's possible.

*/

	is = 1;		ie = Pnx;
	js = Pnx;	je = Pnxny;
	ks = Pnxny;	ke = Pnxnynz;
	if (Pnz == 1) ks = 0, ke = Pnxnynz;

	LOOP3D(I,is,ie,0,Pnxny,0,Pnxnynz,
		Pgrad[I] = (mag[I] > thr) ?
			PWRAP(x[I] - x[I-1]):(yls[I] - yls[I-1]);
	)
	LOOP3D(I,0,Pnx,js,je,0,Pnxnynz,
		Pgrad[I+Pnxnynz] = (mag[I] > thr) ?
			PWRAP(x[I] - x[I-Pnx]):(yls[I] - yls[I-Pnx]);
	)

	/*

	Take divergence (Laplacian of the phase).
	Compute only in the middle of the image (no edge lines/planes).

	*/
	  {	FLOAT	*u,*v,*w;
		u = Pgrad, v = u + Pnxnynz; w = v + Pnxnynz;
		
		is = 1;		ie = Pnx-1;
		js = Pnx;	je = Pnxny - Pnx;
		ks = Pnxny;	ke = Pnxnynz - Pnxny;
		if (Pnz == 1) ks = 0, ke = Pnxnynz;
		
		LOOP3D(I,is,ie,js,je,ks,ke,
		       yls[I] = u[I+1] - u[I];
		)
		LOOP3D(I,is,ie,js,je,ks,ke,
		       yls[I] = yls[I] + v[I+Pnx] - v[I];
	        )
	  }
	
	/* apply Neumann BC for now */

	neubc(yls, Pgrad);		

	cost(yls, Pnx, Pny, Pnz, 1);

	{	register int i,j,k;

		for (k = 0; k < Pnz; k++)
		for (j = 0; j < Pny; j++)
		for (i = 0; i < Pnx; i++) {
			float rhs = -2.*Pndim + Pcosx[i] + Pcosy[j] + Pcosz[k];
			rhs = (i+j+k) ? rhs:1.;		/* 1st eigenvalue = 0 */
			yls[k*Pnxny + j*Pnx + i] /= rhs;
		}
	}

	cost(yls, Pnx, Pny, Pnz, -1);

}

/* 
 * Subroutine	: padjust
 *
 * Description	: adjust the least squares phase to the absolute phase
 *
 */
padjust(y, yls, x, mag, thr)
IMAGE	*y;
FLOAT	*yls;
IMAGE	*x, *mag;
float	thr;
{
	int	length = 16;
	long	*wraps;
	IMAGE	*wspace;			/* work space SMS */
	int	T = (int) (-Pwval), i,loop;
	int	is,ie,js,je,ks,ke;

	/* allocate some spadce */

	if (x == y)				/* if in area = out area */
		wspace = (IMAGE *) (calloc((size_t) (Pnxnynz),
				 (size_t) (sizeof(*wspace))));
	else
		wspace = y;

	wraps = (long *) (calloc((size_t) (length),
				(size_t) (sizeof(*wraps))));

	if (!(wspace && wraps)) {
		fprintf(stderr, "padjust(): cannot allocate arrays\n");
		(void) exit(-1);
	}

#define LOOPADJ(is,ie,js,je,ks,ke,T,I1,I2) {		\
		LOOP3D(I1,is,ie,js,je,ks,ke,		\
			if (( ((float) (abs((int) (T[I1] - T[I2]))	\
			      		)		\
			      )	> Pwval			\
			    ) &&			\
			    ((float) (mag[I1]) > thr	\
			    ) &&			\
			    ((float) (mag[I2]) > thr	\
			    )				\
			   ) wraps[loop] += mag[I1] + mag[I2];	\
		)	\
	}
				/* SMS threshold  */

	is = 1;		ie = Pnx-1;
	js = Pnx;	je = Pnxny - Pnx;
	ks = Pnxny;	ke = Pnxnynz-Pnxny;

	for (loop = 0; loop < length; loop++) {
		for (i = 0; i < Pnxnynz; i++)
			wspace[i] = (IMAGE) (PWRAP(yls[i] - x[i] - T));
		LOOPADJ(is,ie,0,Pnxny,0,Pnxnynz,wspace,I,I-1)
		LOOPADJ(0,Pnx,js,je,0,Pnxnynz,wspace,I,I-Pnx)
		if (Pnz != 1) LOOPADJ(0,Pnx,0,Pnxny,ks,ke,wspace,I,I-Pnxny)
		T = T + Pwvalx2 /(float)(length);
	}

/* find where smallest wraps occur and that's our DC value */

	{
		long minw = wraps[0];
		int  minind = 0;
		float dc_offset;

		for (loop = 0; loop < length; loop++) {
			if (wraps[loop] < minw) {
				minw = wraps[loop];
				minind = loop;
			}
		}

		dc_offset = -Pwval + minind * 2 * Pwval / (float)(length);

		for (i = 0; i < Pnxnynz; i++)
			y[i] = yls[i] -
				(dc_offset + PWRAP(yls[i] - x[i] - dc_offset));

	}


	free(wraps);
	if (x == y) free(wspace);
}

/* 
 * Subroutine	: neubc
 *
 * Description	: takes Laplacian along the boundary assuming Neumann BC
 *
 */
neubc(b, u)
FLOAT	*b, *u;
{
	int	i,j,k,i1,j1,k1,si,sj,sk;
	int	kinc = (Pnz == 1) ? 1:(Pnz-1);
	FLOAT	*v,*w;
	int	is3d = 1;

	v = u + Pnxnynz;
	w = v + Pnxnynz;
	if (Pnz == 1) {
		w = v;		 /* so that code below doesn't blow up */
		is3d = 0;
	}

/*	An example (4x4) of how Neumann BC is appllied.

x =	A	B	C	D
	E	F	G	H
	I	J	K	L
	M	N	O	P

	backward difference			forward difference

dx =	?	B-A	C-B	D-C	dx^2 =	?	OK	OK	?
	?	F-E	G-F	H-G		?	OK	OK	?
	?	J-I	K-J	L-K		?	OK	OK	?
	?	N-M	O-N	P-O		?	OK	OK	?

dy =	?	?	?	?	dy^2 =	?	?	?	?
	E-A	F-B	G-C	H-D		OK	OK	OK	OK
	I-E	J-F	K-G	L-H		OK	OK	OK	OK
	M-I	N-J	O-K	P-L		?	?	?	?

				===>	dx^2	?	?	?	?
					 +   =	?	OK	OK	?
					dy^2	?	OK	OK	?
						?	?	?	?
Where ? (edges), we need:

	(B-A)	0	0	-(D-C)		(E-A)	(F-B)	(G-C)	(H-D)
	(F-E)			-(H-G)	plus	0	0	0	0
	(J-I)			-(L-K)		0	0	0	0
	(N-M)	0	0	-(P-O)		-(M-I)	-(N-J)	-(O-K)	-(P-O)

				plus

	0		(C-B)-(B-A)	(D-C)-(C-B)	0
	(I-E)-(E-A)	0		0		(L-H)-(H-D)
	(M-I)-(I-E)	0		0		(P-L)-(L-H)
	0		(O-N)-(N-M)	(P-O)-(O-N)	0
*/

#define	ELM(X,i,j,k)		X[i + (j)*Pnx + (k)*Pnxny]

#define	FOR_K0			for (k = 0, k1 = 1, sk = 1; k < Pnz;	\
					k += kinc,  k1 = Pnz-1, sk = -1)
#define	FOR_J0			for (j = 0, j1 = 1, sj = 1; j < Pny;	\
					j += Pny-1, j1 = Pny-1, sj = -1)
#define	FOR_I0			for (i = 0, i1 = 1, si = 1; i < Pnx;	\
					i += Pnx-1, i1 = Pnx-1, si = -1)

#define	FOR_K1			for (k = 1; k < Pnz-1; k++)
#define	FOR_J1			for (j = 1; j < Pny-1; j++)
#define	FOR_I1			for (i = 1; i < Pnx-1; i++)

	/* initialize corners */

	FOR_K0	FOR_J0	FOR_I0
		ELM(b,i,j,k) = 	si * ELM(u,i1,j,k) +
				sj * ELM(v,i,j1,k) +
				sk * ELM(w,i,j,k1*is3d) * is3d;
	/* lines */

	FOR_K0	FOR_J0 FOR_I1
		ELM(b,i,j,k) = 	     ELM(u,i+1,j,k) - ELM(u,i,j,k) +
				sj * ELM(v,i,j1,k) +
				sk * ELM(w,i,j,k1*is3d) * is3d;

	FOR_K0	FOR_J1 FOR_I0
		ELM(b,i,j,k) = 	     ELM(v,i,j+1,k) - ELM(v,i,j,k) +
				si * ELM(u,i1,j,k) +
				sk * ELM(w,i,j,k1*is3d) * is3d;

	FOR_K1	FOR_J0	FOR_I0			/* won't exec if 2-D */
		ELM(b,i,j,k) = 	     ELM(w,i,j,k+1) - ELM(w,i,j,k) +
				si * ELM(u,i1,j,k) +
				sj * ELM(v,i,j1,k);

	if (Pnz != 1) {
		/* planes */
		FOR_K0	FOR_J1	FOR_I1
			ELM(b,i,j,k) =	ELM(u,i+1,j,k) - ELM(u,i,j,k) +
					ELM(v,i,j+1,k) - ELM(v,i,j,k) +
					ELM(w,i,j,k1) * sk;
		FOR_K1	FOR_J0	FOR_I1
			ELM(b,i,j,k) =	ELM(u,i+1,j,k) - ELM(u,i,j,k) +
					ELM(w,i,j,k+1) - ELM(w,i,j,k) +
					ELM(v,i,j1,k) * sj;
		FOR_K1	FOR_J1	FOR_I0
			ELM(b,i,j,k) =	ELM(v,i,j+1,k) - ELM(v,i,j,k) +
					ELM(w,i,j,k+1) - ELM(w,i,j,k) +
					ELM(u,i1,j,k) * si;
	}
}
/* 
 * Subroutine	: punwrapinit
 *
 * Description	: initializes the punwrap routine
 *
 */


punwrapinit(wval, nx, ny, nz, niter, thr_type, thrv1, thrv2, mag)
float	wval;
int	nx,ny,nz,niter, thr_type;
float	thrv1, thrv2;
IMAGE	*mag;				/* used if thrv1,2 are computed here */ 
{
	double	pion;
	int	i;

	/* initialize externals */

	Pnx = nx, Pny = ny, Pnz = nz;
	Pnxny = nx*ny; Pnxnynz = nx*ny*nz;
	Pndim = (Pnz == 1) ? 2 : 3;
	Pwval = wval;
	Pwvalx2 = 2. * wval;

	/* allocate */

	Pyls = (FLOAT *) calloc((size_t) (Pnxnynz),
				(size_t) (sizeof(*Pyls)));
	Pgrad = (FLOAT *) calloc((size_t) (Pndim*Pnxnynz),
				(size_t) (sizeof(*Pgrad)));
	Pcosx = (FLOAT *) calloc((size_t) (Pnx),
				(size_t) (sizeof(*Pcosx)));
	Pcosy = (FLOAT *) calloc((size_t) (Pny),
				(size_t) (sizeof(*Pcosy)));
	Pcosz = (FLOAT *) calloc((size_t) (Pnz),
				(size_t) (sizeof(*Pcosz)));
	Pthr = (float *) (calloc((size_t) (niter),
				(size_t) (sizeof(*Pthr))));
	if (!(Pyls && Pgrad && Pcosx && Pcosy && Pcosz && Pthr)) {
		(void) fprintf(stderr,
			"punwrapinit(): cannot allocate arrays\n");
		(void ) exit(-1);
	}

	/* initialize the arrays to compute eigenvalues */

	pion = 4. * atan(1.) / (float) (Pnx);
	for (i = 0; i < Pnx; i++)
		Pcosx[i] = 2. * cos((double) (i*pion));

	pion = 4. * atan(1.) / (float) (Pny);
	for (i = 0; i < Pny; i++)
		Pcosy[i] = 2. * cos((double) (i*pion));

	pion = 4. * atan(1.) / (float) (Pnz);
	for (i = 0; i < Pnz; i++)
		Pcosz[i] = 2. * cos((double) (i*pion));
	if (Pnz == 1) Pcosz[0] = 0.;

	/* initialize the array of thresholds */

	if ((thr_type == 0) || (thr_type == -1)) {
		int maxmag = mag[0];

		for (i = 0; i < Pnxnynz; i++)
			maxmag = maxmag > mag[i] ? maxmag:mag[i];

		if (thr_type == 0) {	/* B_0 map */
			thrv1 = (float) (.05 * maxmag);	/* SMS  5 percent */
			thrv2 = (float) (.15 * maxmag);	/* SMS 15 percent */
							/* 15 works OK */
			printf("punwrap: thresholds %g, %g\n",thrv1,thrv2);
		}
		else {			/* (thr_type == -1) phase contrast */
			thrv1 = (float) (.005 * maxmag);/* SMS .5 percent */
			thrv2 = (float) (.20 * maxmag);	/* SMS 20 percent */
		}
		thr_type = 1;		/* log schedule for both */
	}

	if (niter == 1) Pthr[0] = thrv1;

	else switch(thr_type) {
		case (1):
			for (i = 0; i < niter; i++)	/* log */
				Pthr[i] = thrv1 + (thrv2 - thrv1)
						* log((double) (i+1))
						/ log((double) (niter));
			break;
		case (2):
			for (i = 0; i < niter; i++)	/* linear */
				Pthr[i] = thrv1 + (thrv2 - thrv1)
						* (float) (i)
						/ (float) (niter-1);
			break;
		case (3):
			for (i = 0; i < niter; i++)	/* exp */
				Pthr[i] = thrv1 + (thrv2 - thrv1)
						* exp((double) (i))
						/ exp((double) (niter-1));
			break;
		default:
			fprintf(stdout,"\n");
			fprintf(stdout,
				"punwrapinit(): thr type = %d not valid\n",
				thr_type);
			exit(-1);
			break;
	}
}
/* 
 * Subroutine	: punwrapexit
 *
 * Description	: cleans up the mess
 *
 */
punwrapexit()
{
	free((char *) (Pyls));
	free((char *) (Pgrad));
	free((char *) (Pcosx));
	free((char *) (Pcosy));
	free((char *) (Pcosz));
	free((char *) (Pthr));
	costexit();
}
