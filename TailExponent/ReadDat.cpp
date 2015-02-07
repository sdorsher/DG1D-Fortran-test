#include "ReadDat.h"
using namespace std;


void ReadDat::loadDat(string filename)
{
FILE * fp= fopen(filename,"r");
 double tt;
 double rr;
 double xy;
   while(fscanf(fp,"%lf\t%lf\t%lf\n",&tt,&rr,&xy)==3)
    {
      length++;
      tdat.push_back(tt);
      rdat.push_back(rr);
      xydat.push_back(xy);
    }
fclose(fp);
}
