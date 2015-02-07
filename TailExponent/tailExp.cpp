#include "ReadDat2.h"
#include "calcExponent.h"
#include <iostream>
using namespace std;



int main(void){
  int nmin = 1500;
  double ll=2;
  bool finite = false; // true for r finite, false for xy+
  string filein = "../l2/psi49.dat";
  string fileout = "tailExpl2psi49xy.dat";

  ReadDat2 data=ReadDat2();
  data.loadDat(filein);
  cout << data.tdat[1]<< "\t" << data.rdat[1];
  cout << "\t" << data.xydat[1] << endl;
  int nmax = data.length-10;
  double dt = data.tdat[1]-data.tdat[0];
  double avgalpha;
  if (finite)
    {
      avgalpha =calcExponent(dt, &data.rdat, nmin, nmax, fileout);
    }else
    {
      avgalpha =calcExponent(dt, &data.xydat, nmin, nmax, fileout);
    }

  cout << "The average alpha from t= " << nmin*dt << " to t= " << dt*nmax; 
  cout << " is " << avgalpha << endl;

  if (finite)
    {
      cout << "The theoretical alpha for l = " << ll << " is " << -(2.*ll+3.) << endl;
    }else
    {
      cout << "The theoretical alpha for l = " << ll << " is " << -(ll+2.) << endl;
    }
}


