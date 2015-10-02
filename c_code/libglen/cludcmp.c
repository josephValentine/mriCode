/* Complex form of ludcmp.
   replaces a with LU  decomposition.  indx is output vector of permuatation.
   output d is +/- 1 depending on number of row interchanges being even/odd.
   NOTE: vectors start with index 1, i.e. ar[1] is first element of ar.
*/

#include <math.h>

#define TINY 1.0e-20;

void cludcmp(ar,ai,n,indx,d)
     int n,*indx;
     double **ar,**ai,*d;
{
  int i,imax,j,k;
  double big,dum,sumr,sumi,temp,bigr,bigi,dumr,dumi;
  double *vv,*dvector();
  void nrerror(),free_dvector();
  
  vv=dvector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=hypot(ar[i][j],ai[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sumr=ar[i][j];
      sumi=ai[i][j];
      for (k=1;k<i;k++){
	sumr -= (ar[i][k]*ar[k][j] - ai[i][k]*ai[k][j]);
	sumi -= (ar[i][k]*ai[k][j] + ai[i][k]*ar[k][j]);
      }
      ar[i][j]=sumr;
      ai[i][j]=sumi;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sumr=ar[i][j];
      sumi=ai[i][j];
      for (k=1;k<j;k++) {
	sumr -= (ar[i][k]*ar[k][j] - ai[i][k]*ai[k][j]);
	sumi -= (ar[i][k]*ai[k][j] + ai[i][k]*ar[k][j]);
      }
      ar[i][j]=sumr;
      ai[i][j]=sumi;
      if ( (dum=vv[i]*hypot(sumr,sumi)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=ar[imax][k];
	ar[imax][k]=ar[j][k];
	ar[j][k]=dum;
	dum=ai[imax][k];
	ai[imax][k]=ai[j][k];
	ai[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (hypot(ar[j][j],ai[j][j]) == 0.0) ar[j][j]=TINY;
    if (j != n) { 
      dumr=ar[j][j]/(ar[j][j]*ar[j][j] + ai[j][j]*ai[j][j]);  /* expression for 1 / a[j][j] */
      dumi=(-ai[j][j])/(ar[j][j]*ar[j][j] + ai[j][j]*ai[j][j]);
      for (i=j+1;i<=n;i++) {
	dum = ar[i][j]*dumr - ai[i][j]*dumi;
	ai[i][j] = ar[i][j]*dumi + ai[i][j]*dumr;
	ar[i][j] = dum;
      }
    }
  }
  free_dvector(vv,1,n);
}
  
#undef TINY

