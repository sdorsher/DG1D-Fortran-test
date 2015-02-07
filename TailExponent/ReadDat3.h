#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
class ReadDat3{
 public:
  int length;
  vector<double> tdat,rdat,xydat;
  void loadDat(string filename);
};

