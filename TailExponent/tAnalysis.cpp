#include<math.h>
#include<iostream>
#include<fstream>
using namespace std;
int main(){
  ifstream fp;
  fp.open("../l2/psi25.dat");
  ofstream fp2;
  fp2.open("../l2/tanalysis25.dat");
  float ttmin=0.0;
  float ttmax=0.0;
  float logttmin=0.0;
  float logttmax=0.0;
  float ttlast=0.0;
  float ttdiff=0.0;
  float logttdiff=0.0;
  float logttlast=0.0;
  float logtt=0.0;
  float TOL=1.0*pow(10.0,-5.0);
  int nn=0;
  float tt,psir,psixy;
  char[120] line;
  while(fp.good())
    {
      getline(fp,tt,"\t");
      getline(fp,psir,"\t");
      getline(fp,psixy,"\n");

      if(nn!=0)
	{
	  logtt=log(tt);
	}else
	{
	  logtt=0.0;
	}
      logttdiff=logtt-logttlast;
      ttdiff=tt-ttlast;
      ttmax=max(ttmax,tt);
      logttmax=max(logttmax,logtt);
      ttmin=min(ttmin,tt);
      logttmin=min(logttmin,logtt);
      fp2 << nn<<"\t"<< ttdiff << "\t" << logttdiff << endl;
      nn++;
    }
  fp.close();
  fp2.close();
  cout << ttmax << "," << ttmin <<"\t"<< logttmax << "," << logttmin<<endl;
}
