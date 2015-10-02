#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* 
  resize() takes an x1 by y1 complex matrix with real part m1 and complex part m1i, returns
  in m2 and m2i the fourier interpolation to matrix size x2 by y2.
  If m1i is a null pointer, the input matrix is assumed to be real, and only m2 (not m2i) is
  changed.
*/

resize(m1,m1i,y1,x1,m2,m2i,y2,x2)
double *m1,*m1i,*m2,*m2i;
int x1,y1,x2,y2;

{
  int minx,miny,maxp,evenxin,evenyin,evenxout,evenyout;
  double *tr,*ti,fac,xdist,ydist,xshift,yshift,theta;
  int i,j,org1,org2,real;

  minx = (x1 > x2)?(x2):(x1);
  miny = (y1 > y2)?(y2):(y1);
  maxp = ((x1*y1) > (x2*y2))?(x1*y1):(x2*y2);
  evenxout = !(x2%2);
  evenyout = !(y2%2);
  evenxin = !(x1%2);
  evenyin = !(y1%2);

  if((tr = (double *)malloc(sizeof(double)*maxp)) == NULL) {
    perror("malloc failed in resize");
    exit(1);
  }
  if((ti = (double *)malloc(sizeof(double)*maxp)) == NULL) {
    perror("malloc failed in resize");
    exit(1);
  }
  real = !m1i;
  if(real) {
    if((m1i = (double *)malloc(sizeof(double)*x1*y1)) == NULL) {
      perror("malloc failed in resize");
      exit(1);
    }
    for(i=0;i<x1*y1;i++)
      m1i[i] = 0.0;
    if((m2i = (double *)malloc(sizeof(double)*x2*y2)) == NULL) {
      perror("malloc failed in resize");
      exit(1);
    }
  }
  fft2d(m1,m1i,tr,ti,y1,x1,0);
  oddfftshift2d(tr,y1,x1);
  oddfftshift2d(ti,y1,x1);
  for(i=0;i<x2*y2;i++) {
    m2[i] = 0.0;
    m2i[i] = 0.0;
  }

  xshift = -(double)x1 / (double)x2 / 2.0;
  yshift = (double)y1 / (double)y2 / 2.0;

  xshift += 0.5;
  yshift -= 0.5;

  org1 = x1/2 - minx/2 + x1 * (y1/2 - miny/2);
  org2 = x2/2 - minx/2 + x2 * (y2/2 - miny/2);

  for(j=0; j<miny; j++) { 
    ydist = (double)((miny / 2) - j)/(double)y1;
    for(i=0; i<minx; i++) {
      xdist = (double)(i - (minx / 2))/(double)x1;
      theta = 2.0 * M_PI * (xshift * xdist + yshift * ydist);

      m2[j*x2+i+org2] = tr[j*x1+i+org1] * cos(theta) - ti[j*x1+i+org1] * sin(theta);
      m2i[j*x2+i+org2] = ti[j*x1+i+org1] * cos(theta) + tr[j*x1+i+org1] * sin(theta);
      if(!j) {
	if(evenyin && (y2 > y1)) {
	  m2[i+org2] *= 0.5;
	  m2i[i+org2] *= 0.5;
	  theta = 2.0 * M_PI * (xshift * xdist - yshift * ydist);
	  m2[i+org2+miny*x2] = 0.5 * (tr[i+org1] * cos(theta) - ti[i+org1] * sin(theta));
	  m2i[i+org2+miny*x2] = 0.5 * (ti[i+org1] * cos(theta) + tr[i+org1] * sin(theta));
	}
	else if(evenyout && (y2 < y1)) {
	  theta = 2.0 * M_PI * (xshift * xdist - yshift * ydist);
	  m2[i+org2] += tr[i+org1+miny*x1] * cos(theta) - ti[i+org1+miny*x1] * sin(theta);
	  m2i[i+org2] += ti[i+org1+miny*x1] * cos(theta) + tr[i+org1+miny*x1] * sin(theta);
	}
      }
      if(!i) {
	if(evenxin && (x2 > x1)) {
	  m2[j*x2+org2] *= 0.5;
	  m2i[j*x2+org2] *= 0.5;
	  theta = 2.0 * M_PI * (-xshift * xdist + yshift * ydist);
	  m2[j*x2+org2+minx] = 0.5 * (tr[j*x1+org1] * cos(theta) - ti[j*x1+org1] * sin(theta));
	  m2i[j*x2+org2+minx] = 0.5 * (ti[j*x1+org1] * cos(theta) + tr[j*x1+org1] * sin(theta));
	}
	else if(evenxout && (x2 < x1)) {
	  theta = 2.0 * M_PI * (-xshift * xdist + yshift * ydist);
	  m2[j*x2+org2] += tr[j*x1+org1+minx] * cos(theta) - ti[j*x1+org1+minx] * sin(theta);
	  m2i[j*x2+org2] += ti[j*x1+org1+minx] * cos(theta) + tr[j*x1+org1+minx] * sin(theta);
	}
      }
    }
    
  }
/* four points may need to be fixed */

  xdist = (double)(x2 / 2)/(double)x1;
  ydist = (double)(y2 / 2)/(double)y1;
  theta = 2.0 * M_PI * (xshift * xdist - yshift * ydist);
  if(x1 > x2) {
    if(y1 > y2) {
      if(evenxout && evenyout) {
	m2[org2] += tr[org1 + miny*x1 + minx] * cos(theta) - ti[org1 + miny*x1 + minx] * sin(theta);
	m2i[org2] += ti[org1 + miny*x1 + minx] * cos(theta) + tr[org1 + miny*x1 + minx] * sin(theta);
      }
    }
    if(y1 < y2) {
      if(evenxout && evenyin) {
	m2[org2 + miny*x2] += 0.5 * (tr[org1+minx] * cos(theta) - ti[org1 + minx] * sin(theta));
	m2i[org2 + miny*x2] += 0.5 * (ti[org1+minx]* cos(theta) + tr[org1 + minx] * sin(theta));
      }
    }
  }
  if(x1 < x2) {
    if(y1 > y2) {
      if(evenyout && evenxin) {
	m2[org2 + minx] += 0.5 * (tr[org1+miny*x1] * cos(theta) - ti[org1+miny*x1] * sin(theta));
	m2i[org2 + minx] += 0.5 * (ti[org1+miny*x1] * cos(theta) + tr[org1+miny*x1] * sin(theta));
      }
    }
    if(y1 < y2) {
      if(evenxin && evenyin) {
	m2[org2 + miny*x2] *= 0.5;
	m2i[org2 + miny*x2] *= 0.5;
	m2[org2+miny*x2+minx] = 0.25 * (tr[org1]*cos(theta) - ti[org1]*sin(theta));
	m2i[org2+miny*x2+minx] = 0.25 * (ti[org1]*cos(theta) + tr[org1]*sin(theta));
      }
    }
  }

  fftshift2d(m2,y2,x2);
  fftshift2d(m2i,y2,x2);
  fft2d(m2,m2i,m2,m2i,y2,x2,1);
  fac = (double)(x2 * y2) / (double)(x1 * y1);
  for(i=0; i<x2*y2; i++) {
    m2[i] *= fac;
    m2i[i] *= fac;
  }
  if(real) {
    free(m1i);
    free(m2i);
  }
  free(tr);
  free(ti);
}



