#include "calcExponent.h"
double calcExponent(double dt, vector<double>* xx, int nmin, int nmax, string outfile)
{
  double thisalpha,sumalpha=0., avgalpha;
  vector<double> dxdt = derivCentral4thOrder(dt, xx);
  

  ofstream fs;
  fs.open(outfile);
  for(int i=nmin; i<=nmax; i++)
    {
      thisalpha = i*dt*dxdt[i]/(*xx)[i];
      sumalpha+=thisalpha;
      fs << i*dt << "\t" << thisalpha << endl;
    }
  avgalpha=sumalpha/(nmax-nmin+1.0);
  fs.close();

  return avgalpha;
}

