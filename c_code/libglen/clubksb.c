/* Complex version.
   solves AX=b.  a in LU decomposition. indx is permutation vector.
   solution overwrites b.
   NOTE: vectors start with index 1.
*/

#include<math.h>

void clubksb(ar,ai,n,indx,br,bi)
double **ar,**ai,br[],bi[];
int n,*indx;
{
  int i,ii=0,ip,j;
  double sumr,sumi;
  double dum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sumr=br[ip];
    sumi=bi[ip];
    br[ip]=br[i];
    bi[ip]=bi[i];
    if (ii)
      for (j=ii;j<=i-1;j++){
	sumr -= (ar[i][j]*br[j] - ai[i][j]*bi[j]);
	sumi -= (ar[i][j]*bi[j] + ai[i][j]*br[j]);
      }
    else if (hypot(sumr,sumi)) ii=i;
    br[i]=sumr;
    bi[i]=sumi;
  }
  for (i=n;i>=1;i--) {
    sumr=br[i];
    sumi=bi[i];
    for (j=i+1;j<=n;j++) {
      sumr -= (ar[i][j]*br[j] - ai[i][j]*bi[j]);
      sumi -= (ar[i][j]*bi[j] + ai[i][j]*br[j]);
    }
    br[i] = (sumr*ar[i][i] + sumi*ai[i][i]) / (ar[i][i]*ar[i][i] + ai[i][i]*ai[i][i]);
    bi[i] = (sumi*ar[i][i] - sumr*ai[i][i]) / (ar[i][i]*ar[i][i] + ai[i][i]*ai[i][i]);
  }
}

