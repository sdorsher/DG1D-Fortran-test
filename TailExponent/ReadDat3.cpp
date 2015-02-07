#include "ReadDat3.h"
using namespace std;
//compile with -std=c++11

void ReadDat3::loadDat(string filename)
{
  ifstream fp;
  fp.open(filename);
  double data;
  int num=0;
  while(fp>>data)
    {
      switch (num%3)
	{
	case 0: 
	  tdat.push_back(data);
	  break;
	case 1:
	  rdat.push_back(data);
	  break;
	case 2:
	  xydat.push_back(data);
	  length++;
	  break;
	}
      num++;
    }
  fp.close();
}
