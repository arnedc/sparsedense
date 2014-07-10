#include <stdio.h>
#include <stdlib.h>
//#include <mkl_lapack.h>
#include "shared_var.h"
//#include <mkl_blas.h>
#include "CSRdouble.hpp"
#include "ParDiSO.hpp"
#include "RealMath.hpp"
#include <cassert>

extern "C" {
    void dgemm_ ( const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc );
    void dpotrf_( const char* uplo, const int* n, double* a, const int* lda, int* info );
    void dpotrs_( const char* uplo, const int* n, const int* nrhs,const double* a, const int* lda, double* b, const int* ldb, int* info );
}

/**
 * @brief I tried to create the Schur complement only using sparse B's and the Schur complement modus of PARDISO, but it did not work due to the fact that C can get rectangular
 *
 * @param A Full sparse (0,0)-block of which we want the Schur complement in matrix C
 * @param BT_i Sparse (1,0)-block of C corresponding to T_ij
 * @param B_j Sparse (0,1)-block of C corresponding to T_ij
 * @param T_ij Dense (1,1)-block of C
 * @param lld_T local leading dimension of T_ij
 * @return int
 **/
int make_Sij_sparse_parallel(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T) {
    CSRdouble C,S;
    if (iam==0) {
        cout << "***                                           [ A      B_j ] *** " << endl;
        cout << "***                                           [            ] *** " << endl;
        cout << "*** G e n e r a t i n g    m a t r i x    C = [            ] *** " << endl;
        cout << "***                                           [    t       ] *** " << endl;
        cout << "***                                           [ B_i    T_ij] *** " << endl;
    }
    create2x2BlockMatrix_denseT_lldT(A, BT_i, B_j, T_ij, lld_T,  C);
    blacs_barrier_ ( &ICTXT2D, "A" );

    S.nrows=BT_i.nrows;
    S.ncols=B_j.ncols;

    calculateSchurComplement(C, 11, S);
    blacs_barrier_ ( &ICTXT2D, "A" );

    CSR2dense_lld ( S, T_ij, lld_T ) ;
    return 0;
}

/**
 * @brief BT_i and B_j are converted to dense matrices in each process to solve the sparse system AX=B_j and afterwards do BT_i * X. X is stored as a dense matrix in AB_sol
 *
 * @param A Sparse (0,0)-block of which we want to compute the Schur complement in matrix C
 * @param BT_i Sparse (1,0)-block of C corresponding to T_ij.
 * @param B_j Sparse (0,1)-block of C corresponding to T_ij
 * @param T_ij Dense (1,1)-block of C
 * @param lld_T local leading dimension of T_ij
 * @param AB_sol_out Dense solution of AX=B_j (output)
 * @return int
 **/
int make_Sij_parallel_denseB(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T, double * AB_sol_out) {

    double *BT_i_dense, *B_j_dense, *AB_sol;

    assert(A.nrows == BT_i.ncols);

    BT_i_dense=(double *) calloc(BT_i.nrows * BT_i.ncols,sizeof(double));
    B_j_dense=(double *) calloc(B_j.nrows * B_j.ncols,sizeof(double));
    AB_sol=(double *) calloc(A.nrows * B_j.ncols,sizeof(double));

    CSR2dense(BT_i,BT_i_dense);
    CSR2dense(B_j,B_j_dense);

    solveSystem(A, AB_sol_out,B_j_dense, 2, B_j.ncols);

    printf("Processor %d finished solving system AX=B\n",iam);

    dgemm_("N","N",&(BT_i.nrows),&(B_j.ncols),&(BT_i.ncols),&d_negone,BT_i_dense,&(BT_i.nrows),
           AB_sol_out,&(A.nrows),&d_one,T_ij,&lld_T
          );

    return 0;
}
