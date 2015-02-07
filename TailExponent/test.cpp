#include "ReadDat3.h"
#include <iostream>
using namespace std;



int main(void){
  ReadDat3 l2psi25dat=ReadDat3();
  l2psi25dat.loadDat("testdat.dat");
  cout << l2psi25dat.tdat[1]<< "\t" << l2psi25dat.rdat[1];
  cout << "\t" << l2psi25dat.xydat[1] << endl;
}
