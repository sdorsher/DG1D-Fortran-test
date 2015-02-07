#include <iostream>
#include <math.h>
#include <vector>
#include <string>
using namespace std;
class ReadDat{
 public:
  int length;
  vector<double> tdat,rdat,xydat;
  void loadDat(string filename);
};

