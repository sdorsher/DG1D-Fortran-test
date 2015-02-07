#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

int main(){
vector<float> tdat,rdat,xydat;
FILE * fp= fopen("testdat.dat","r");
 float tt;
 float rr;
 float xy;
   while(fscanf(fp,"%f\t%f\t%f\n",&tt,&rr,&xy)==3)
    {
      tdat.push_back(tt);
      rdat.push_back(rr);
      xydat.push_back(xy);
    }
fclose(fp);
cout << tdat[1] << "\t" << rdat[1] << "\t" << xydat[1] << endl;
}
