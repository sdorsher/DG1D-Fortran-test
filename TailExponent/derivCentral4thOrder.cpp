#include "derivCentral4thOrder.h"
using namespace std;

vector<double>  derivCentral4thOrder(double dt, vector<double>* xx)
{
  int nmax = xx->size();
  vector<double> dxbydt(nmax, 0.0);
  double coeff[5]={1.0/12.0, -2.0/3.0, 0.0, 2./3., -1./12.};
  //handle boundaries
  int order = 4;
  for(int i = 0; i<nmax; i++) 
    { 
      for(int k=0; k<=order; k++)
	{
	  if((k+i>=order/2)&&(k+i<nmax+order/2))
	    {
	      dxbydt[i]+=coeff[k]*(*xx)[i-order/2+k]/dt;
	    }
	}
    }
  return dxbydt;
}

