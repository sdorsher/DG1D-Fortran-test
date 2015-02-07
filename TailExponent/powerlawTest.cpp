#include "calcExponent.h"
#include <iostream>
#include <math.h>

using namespace std;


int main(void) {


  vector<double> xx(10000,0.0);
  double dt = 0.01;
  int nmin = 1500;
  int nmax = 10000-10;
  double ll=2;

  for(int i=0; i< 10000; i++)
    {
      xx[i]=pow((dt*i),-(2.*ll+3.));
    }
	

  double avgalpha =calcExponent(dt, &xx, nmin, nmax, "testPowerLaw.dat");


  cout << "The average alpha from t= " << nmin*dt << " to t= " << dt*nmax; 
  cout << " is " << avgalpha << endl;
  cout << "The theoretical alpha for l = " << ll << " is " << -(2.*ll+3.) << endl;
}
