#ifndef config_hpp
#define config_hpp

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>

using std::complex;
using std::endl;
using std::setw;
using std::cout;
using std::fstream;
using std::ios;
using std::string;
using std::vector;

enum MatrixStorage
{
  NOT_SET   = 0,
  NORMAL    = 1,
  TRANSPOSE = 2,
  SYMMETRIC = 3,  
  HERMITIAN = 4,  
};

#endif
