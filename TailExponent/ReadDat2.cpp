#include "ReadDat2.h"
using namespace std;


void ReadDat2::loadDat(string filename)
{
  ifstream fp;
  fp.open(filename);
  double tt;
  double rr;
  double xy;
  string line;
  while(fp.good())
    {
      getline(fp,line);
      sscanf(line.c_str(),"%lf\t%lf\t%lf\n",&tt,&rr,&xy);
      length++;
      tdat.push_back(tt);
      rdat.push_back(rr);
      xydat.push_back(xy);
    }
  fp.close();
}
