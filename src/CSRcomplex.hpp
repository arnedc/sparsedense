#ifndef CSRcomplex_hpp
#define CSRcomplex_hpp

// include files
// =============
#include "config.hpp"
// =============

class CSRcomplex
{
public:
    string            name;
    int               nrows;
    int               ncols;
    int               nonzeros;
    int*              pRows;
    int*              pCols;
    complex<double>*  pData;

    MatrixStorage     matrixType;

public:
    CSRcomplex();
    ~CSRcomplex();

    void  allocate(int n, int m, int nzeros);
    void  residual(complex<double>* r, complex<double>* x, complex<double>* b);
    void  multiply(complex<double>* x, complex<double>* y);
    void  multiplyH(complex<double>* x, complex<double>* b);
    void  multiplyS(complex<double>* x, complex<double>* y);
    void  multiplyN(complex<double>* x, complex<double>* y);
    void  multiplyT(complex<double>* x, complex<double>* y);
    void  sortColumns();
    void  make(int n, int m, int nzeros, int* prows, int* pcols, complex<double>* pdata);
    void  writeToFile(const char* filename, ios::openmode mode = ios::out) const;
    void  loadFromFile(const char* file, ios::openmode mode = ios::out);
    void  savedebug(const char* filename) const;
    void  clear();
};


#endif
