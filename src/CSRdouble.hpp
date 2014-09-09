#ifndef CSRdouble_hpp
#define CSRdouble_hpp

// include files
// =============
#include "config.hpp"
// =============

class CSRdouble {
public:
    string         name;
    int            nrows;
    int            ncols;
    int            nonzeros;
    int*           pRows;
    int*           pCols;
    double*        pData;
    MatrixStorage  matrixType;

public:
    CSRdouble();
    ~CSRdouble();

    void  allocate ( int n, int m, int nzeros );
    void  loadFromFile ( const char* file, ios::openmode mode = ios::out );
    void  loadFromFileCOO ( const char* file );
    void  make ( int n, int m, int nzeros, int* prows, int* pcols, double* pdata );
    void  make2 ( int n, int m, int nzeros, int* prows, int* pcols, double* pdata );
    void  transposeIt ( int block_size );

    void  addBCSR ( CSRdouble& B );
    void  extendrows ( CSRdouble& B, int startrowB, int nrowsB );

    void  residual ( double* r, double* x, double* b );
    void  multiply ( double* x, double* y );
    void  multiplyS ( double* x, double* y );
    void  multiplyN ( double* x, double* y );
    void  multiplyT ( double* x, double* y );
    void  sortColumns();
    void  fillSymmetric();
    void  reduceSymmetric();
    void  writeToFile ( const char* filename, ios::openmode mode = ios::out ) const;
    void  savedebug ( const char* filename ) const;
    void  clear();
};


#endif
