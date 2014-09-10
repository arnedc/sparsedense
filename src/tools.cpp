#include "CSRdouble.hpp"
#include "config.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <cstring>
#include "shared_var.h"

void printdense ( int m, int n, double *mat, char *filename ) {
    FILE *fd;
    fd = fopen ( filename,"w" );
    if ( fd==NULL )
        printf ( "error creating file" );
    int i,j;
    for ( i=0; i<m; ++i ) {
        //fprintf ( fd,"[\t" );
        for ( j=0; j<n; ++j ) {
            fprintf ( fd,"%12.8g\t",* ( mat+i*n +j ) );
        }
        fprintf ( fd,"\n" );
    }
    fclose ( fd );
}

//converting a dense matrix (m x n) stored column-wise to CSR format
void dense2CSR ( double *mat, int m, int n, CSRdouble& A ) {
    int i,j, nnz;
    double *pdata;
    int  *prows,*pcols;

    nnz=0;

    for ( i=0; i<m; ++i ) {
        for ( j=0; j<n; ++j ) {
            if ( abs ( * ( mat+i*n+j ) ) >1e-10 ) {
                nnz++;
            }
        }
    }

    prows= new int [m+1];
    if ( prows == NULL ) {
        printf ( "unable to allocate memory for prows in dense2CSR (required: %ld bytes)\n", (m+1) * sizeof ( int ) );
        exit(1);
    }
    pcols= new int [nnz];
    if ( pcols == NULL ) {
        printf ( "unable to allocate memory for pcols in dense2CSR (required: %ld bytes)\n", nnz * sizeof ( int ) );
        exit(1);
    }
    pdata= new double [nnz];
    if ( pdata == NULL ) {
        printf ( "unable to allocate memory for pdata in dense2CSR (required: %ld bytes)\n", nnz * sizeof ( double ) );
        exit(1);
    }

    *prows=0;
    nnz=0;
    for ( i=0; i<m; ++i ) {
        for ( j=0; j<n; ++j ) {
            if ( abs ( * ( mat+j*m+i ) ) >1e-10 ) { //If stored column-wise (BLAS), then moving through a row is going up by m (number of rows).
                * ( pdata+nnz ) =* ( mat+j*m+i );
                * ( pcols+nnz ) =j;
                nnz++;
            }
        }
        * ( prows+i+1 ) =nnz;
    }
    A.clear();
    A.make ( m,n,nnz,prows,pcols,pdata );
}

//converting a dense matrix (m x n) stored column-wise to CSR format at specific submatrix
void dense2CSR_sub ( double *mat, int m, int n, int lld_mat, CSRdouble& A, int startrow, int startcol ) {
    int i,j, nnz, rows, cols;
    double *pdata;
    int  *prows,*pcols;

    assert(A.nrows>=startrow + m);
    assert(A.ncols>=startcol + n);

    nnz=0;

    for ( i=0; i<m; ++i ) {
        for ( j=0; j<n; ++j ) {
            if ( abs ( * ( mat+j*lld_mat+i ) ) >1e-10 ) {
                nnz++;
            }
        }
    }

    prows= new int [A.nrows + 1];
    if ( prows == NULL ) {
        printf ( "unable to allocate memory for prows in dense2CSR (required: %ld bytes)\n", (A.nrows+1) * sizeof ( int ) );
        exit(1);
    }
    pcols= new int[nnz];
    if ( pcols == NULL ) {
        printf ( "unable to allocate memory for pcols in dense2CSR (required: %ld bytes)\n", nnz * sizeof ( int ) );
        exit(1);
    }
    pdata= new double[nnz];
    if ( pdata == NULL ) {
        printf ( "unable to allocate memory for pdata in dense2CSR (required: %ld bytes)\n", nnz * sizeof ( double ) );
        exit(1);
    }

    *prows=0;
    nnz=0;
    for (i=1; i<=startrow; i++) {
        * (prows+i)=0;
    }
    for ( i=0; i<m; ++i ) {
        for ( j=0; j<n; ++j ) {
            if ( abs ( * ( mat+j*lld_mat+i ) ) >1e-10 ) { //If stored column-wise (BLAS), then moving through a row is going up by lld_mat (number of rows).
                * ( pdata+nnz ) = * ( mat+j*lld_mat+i );
                * ( pcols+nnz ) = j+startcol;
                nnz++;
            }
        }
        * ( prows+i+startrow+1 ) =nnz;
    }
    for (i=startrow+m+1; i<=A.nrows; ++i) {
        *(prows+i)=nnz;
    }
    rows=A.nrows;
    cols=A.ncols;
    A.clear();
    A.make ( rows,cols,nnz,prows,pcols,pdata );
}

//Convert CSR to column-wise stored dense matrix
void CSR2dense ( CSRdouble& matrix,double *dense ) {
    int i, row;
    row=0;
    for ( i=0; i<matrix.nonzeros; ++i ) {
        while ( i==matrix.pRows[row+1] )
            row++;
        * ( dense + row + matrix.nrows * matrix.pCols[i] ) =matrix.pData[i] ;
    }
}

void CSR2dense_lld ( CSRdouble& matrix,double *dense, int lld_dense ) {
    int i, row;
    row=0;
    for ( i=0; i<matrix.nonzeros; ++i ) {
        while ( i==matrix.pRows[row+1] )
            row++;
        * ( dense + row + lld_dense * matrix.pCols[i] ) =matrix.pData[i] ;
    }
}

void mult_CSRA_denseB_storeCSRC ( CSRdouble& A, double *B, bool trans,
                                  CSRdouble& C ) {
    int i, j,row, col, C_nnz, C_ncols, *prows;
    double cij;

    C_ncols=C.ncols;

    vector<int> Ccols;
    vector<double> Cdata;

    prows = new int[A.nrows + 1];
    C_nnz = 0;

    for ( row=0; row<A.nrows; ++row ) {
        for ( col=0; col<C.ncols; ++col ) {
            cij=0;
            for ( i=A.pRows[row]; i<A.pRows[row+1]; ++i ) {
                j = A.pCols[i];
                cij += A.pData[i] * * ( B + j + A.ncols * col ) ;
            }
            if ( abs ( cij ) >1e-10 ) {
                C_nnz++;
                Ccols.push_back ( col );
                Cdata.push_back ( cij );
            }
        }
        prows[row+1]=C_nnz;
    }
    double* pdata = new double[C_nnz];
    int*    pcols = new int[C_nnz];

    memcpy ( pcols, &Ccols[0], C_nnz*sizeof ( int ) );
    memcpy ( pdata, &Cdata[0], C_nnz*sizeof ( double ) );

    Ccols.clear();
    Cdata.clear();
    
    C.clear();

    C.make ( A.nrows,C_ncols,C_nnz,prows,pcols,pdata );

    if ( trans )
        C.transposeIt ( 1 );
}

/**
 * @brief Computes some columns of sparse A with dense B and stores it into sparse C starting at a given column. A and C must have the same number of rows
 *
 * @param A Sparse matrix of which some columns are selected to be multiplied with B
 * @param B Dense matrix to be multiplied with the columns of B (stored row-wise)
 * @param lld_B local leading dimension of B (should always be larger than Cncols)
 * @param Acolstart First column of the submatrix of A to be multiplied with B (if Acolstart > A.ncols no multiplication is performed, but no error is returned)
 * @param Ancols Number of columns in the submatrix of A to be multiplied with B (if Acolstart + Ancols > A.ncols no multiplication is performed for columns > a.ncols, but no error is returned
 * @param Ccolstart First column of the submatrix of C where the result of A * B is stored.
 * @param Cncols Number of columns of the submatrix of C which are calculated in A * B.
 * @param C Sparse matrix containing the result of A * B (output)
 * @param trans Can be set to 1 of we want to store the transposed result B' * A'.
 * @return void
 **/
void mult_colsA_colsC ( CSRdouble& A, double *B, int lld_B, int Acolstart, int Ancols, int Ccolstart, int Cncols, //input
                        CSRdouble& C, bool trans ) {
    int i, j,row, col, C_nnz,C_ncols, *prows;
    double cij;

    /*assert(Cncols < lld_B);
    assert(Ccolstart+Cncols <= C.ncols);*/

    C_ncols=C.ncols;

    vector<int> Ccols;
    vector<double> Cdata;

    prows = new int[A.nrows + 1];
    C_nnz = 0;
    prows[0]=0;

    for ( row=0; row<A.nrows; ++row ) {
        for ( col=Ccolstart; col<Ccolstart+Cncols; ++col ) {
            cij=0;
            for ( i=A.pRows[row]; i<A.pRows[row+1]; ++i ) {
                j = A.pCols[i];
                if ( j>=Acolstart && j<Acolstart+Ancols )
                    cij += A.pData[i] * * ( B + col-Ccolstart + lld_B * ( j-Acolstart ) ) ;
            }
            if ( abs ( cij ) >1e-10 ) {
                C_nnz++;
                Ccols.push_back ( col );
                Cdata.push_back ( cij );
            }
        }
        prows[row+1]=C_nnz;
    }
    double* pdata = new double[C_nnz];
    int*    pcols = new int[C_nnz];

    memcpy ( pcols, &Ccols[0], C_nnz*sizeof ( int ) );
    memcpy ( pdata, &Cdata[0], C_nnz*sizeof ( double ) );

    Ccols.clear();
    Cdata.clear();
    C.clear();

    C.make ( A.nrows,C_ncols,C_nnz,prows,pcols,pdata );

    if ( trans )
        C.transposeIt ( 1 );
}

void mult_colsA_colsC_denseC ( CSRdouble& A,double *B, int lld_B, int Acolstart, int Ancols, int Ccolstart, int Cncols,
                               double *C, int lld_C, bool sum, double alpha ) {
    int i, j,row, col;

    /*assert(Cncols < lld_B);
    assert(Ccolstart+Cncols <= C.ncols);*/

    for ( row=0; row<A.nrows; ++row ) {
        for ( col=Ccolstart; col<Ccolstart+Cncols; ++col ) {
            if (!sum)
                *(C + row + col * lld_C) = 0;
            for ( i=A.pRows[row]; i<A.pRows[row+1]; ++i ) {
                j = A.pCols[i];
                if ( j>=Acolstart && j<Acolstart+Ancols )
                    *(C + row + col * lld_C) = *(C + row + col * lld_C) + alpha * A.pData[i] * * ( B + lld_B * (col-Ccolstart) + j-Acolstart  ) ;
            }
        }
    }
}

/**
 * @brief This function adds sparse matrix B to the CSRdouble sparse matrix A
 *
 * @param B Sparse matrix to be added to the CSRdouble matrix
 * @return void
 **/
void CSRdouble::addBCSR ( CSRdouble& B ) {
    double sum;
    int Acolindex, Bcolindex,nonzeroes, colindex;
    int * ABprows = new int[nrows+1];
    ABprows[0]=0;
    nonzeroes=0;

    vector<int> ABcols;
    vector<double> ABdata;

    if ( nrows != B.nrows ) {
        printf ( "rows of A (%d) are not the same as rows of B (%d)",nrows,B.nrows );
    }
    if ( ncols != B.ncols ) {
        printf ( "rows of A (%d) are not the same as rows of B (%d)",ncols,B.ncols );
    }
    assert ( nrows == B.nrows );
    assert ( ncols == B.ncols );

    for ( int i=0; i<nrows; ++i ) {
        Acolindex= pRows[i];
        Bcolindex=B.pRows[i];

        while ( Acolindex < pRows[i+1] && Bcolindex < B.pRows[i+1] ) {
            if ( pCols[Acolindex] == B.pCols[Bcolindex] ) {
                colindex=pCols[Acolindex];
                sum=pData[Acolindex] + B.pData[Bcolindex];
                ABdata.push_back ( sum );
                ABcols.push_back ( colindex );
                Acolindex++, Bcolindex++, nonzeroes++;
            } else if ( pCols[Acolindex] < B.pCols[Bcolindex] ) {
                colindex=pCols[Acolindex];
                sum=pData[Acolindex];
                ABdata.push_back ( sum );
                ABcols.push_back ( colindex );
                Acolindex++, nonzeroes++;
            } else {
                colindex=B.pCols[Bcolindex];
                sum=B.pData[Bcolindex];
                ABdata.push_back ( sum );
                ABcols.push_back ( colindex );
                Bcolindex++, nonzeroes++;
            }
        }
        while ( Acolindex < pRows[i+1] ) {
            colindex=pCols[Acolindex];
            sum=pData[Acolindex];
            ABdata.push_back ( sum );
            ABcols.push_back ( colindex );
            Acolindex++, nonzeroes++;
        }
        while ( Bcolindex < B.pRows[i+1] ) {
            colindex=B.pCols[Bcolindex];
            sum=B.pData[Bcolindex];
            ABdata.push_back ( sum );
            ABcols.push_back ( colindex );
            Bcolindex++, nonzeroes++;
        }
        ABprows[i+1]=nonzeroes;
    }
    double* pdata = new double[nonzeroes];
    int*    pcols = new int[nonzeroes];

    if ( ABprows[nrows]!=nonzeroes )
        printf ( "last element of prows (%d) not equal to number of nonzeroes (%d)\n",ABprows[nrows],nonzeroes );

    memcpy ( pcols, &ABcols[0], nonzeroes*sizeof ( int ) );
    memcpy ( pdata, &ABdata[0], nonzeroes*sizeof ( double ) );
    ABcols.clear();
    ABdata.clear();

    delete[] pRows;
    delete[] pCols;
    delete[] pData;

    make ( B.nrows, B.ncols, nonzeroes, ABprows, pcols, pdata );
}

/**
 * @brief Extends the sparse CSRDouble with a certain number of rows from sparse matrix B. The rows are simply added to the end of the original matrix
 *
 * @param B Sparse matrix contianing the rows to be added to the original matrix, must have identical number of colums as original matrix
 * @param startrowB First row of B to be added to the original (must be smaller than number of rows in B
 * @param nrowsB Number of rows to be added to the original (if nrowsB + startrowB exceeds B.nrows, empty rows are added to the original matrix)
 * @return void
 **/
void CSRdouble::extendrows ( CSRdouble& B, int startrowB, int nrowsB ) {
    assert ( ncols==B.ncols );
    assert (startrowB < B.nrows);
    int  nonzeroes, nonzeroesB, i, colindex;
    int  n        = nrows + nrowsB;
    int* prows;
    int* pcols;
    double* pdata;

    if (startrowB+nrowsB > B.nrows) {
        nonzeroesB = B.nonzeros - B.pRows[startrowB];
        n = B.nrows - startrowB + nrows;
    }
    else
        nonzeroesB = B.pRows[startrowB + nrowsB] - B.pRows[startrowB];
    nonzeroes=nonzeroesB + nonzeros;
    colindex=B.pRows[startrowB];

    prows = new int [n+1];
    pcols = new int [nonzeroes];
    pdata = new double [nonzeroes];

    memcpy ( prows, & (pRows[0] ), nrows * sizeof(int));
    memcpy ( pcols, & (pCols[0] ), nonzeros * sizeof(int));
    memcpy ( pdata, & (pData[0] ), nonzeros * sizeof(double));
    for (i=0; i<=nrowsB; ++i) {
        if(startrowB+i <= B.nrows)
            prows[nrows+i]= nonzeros + B.pRows[startrowB+i]-B.pRows[startrowB];
    }

    assert(prows[n]==nonzeroes);

    if(prows[n] != nonzeroes)
        printf("nonzeroes (%d) not equal to last element of prows (%d)", nonzeroes, prows[n]);

    memcpy ( &(pcols[nonzeros]), & ( B.pCols[colindex] ), nonzeroesB * sizeof(int));
    memcpy ( &(pdata[nonzeros]), & ( B.pData[colindex] ), nonzeroesB * sizeof(double));

    delete[] pRows;
    delete[] pCols;
    delete[] pData;

    make ( n, B.ncols, nonzeroes, prows, pcols, pdata );
}
