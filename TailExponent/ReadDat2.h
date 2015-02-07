#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
class ReadDat2{
 public:
  int length;
  vector<double> tdat,rdat,xydat;
  void loadDat(string filename);
};

