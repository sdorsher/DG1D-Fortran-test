#include "ReadDat.h"
#include "derivCentral4thOrder.h"
#include <iostream>
#include <fstream>
using namespace std;



int main(void){
  double dt = 0.01;
  double slope = -3.0;
  vector<double>  tt(1000, 0.0),linear(1000, 0.0);
  for(int i=0; i<1000; i++)
    {
      tt[i] = dt*i;
      linear[i]=slope*tt[i];
    }

  
  vector<double> dxdt = derivCentral4thOrder(dt, &linear);
  ofstream fs;
  fs.open("testDeriv.dat");
  for(int i=0; i<1000; i++)
    {
      fs << tt[i] << "\t" << linear[i] << "\t" << dxdt[i] << endl;
    }

  fs.close();
  
}
