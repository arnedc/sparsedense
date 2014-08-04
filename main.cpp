#include <stdio.h>
#include <stdlib.h>
#include "src/shared_var.h"
#include <mkl_blas.h>
#include <string.h>
#include <mkl_lapack.h>
#include <mkl_scalapack.h>
#include <mkl_pblas.h>
#include "src/config.hpp"
#include "src/CSRdouble.hpp"
#include "src/IO.hpp"
#include "src/ParDiSO.hpp"
#include "src/RealMath.hpp"
#include "src/timing.hpp"
//#include <mpi.h>
#include "src/smat.h"


extern "C" {
    int MPI_Init(int *, char ***);
    int MPI_Dims_create(int, int, int *);
    void blacs_pinfo_ ( int *mypnum, int *nprocs );
    void blacs_setup_ ( int *mypnum, int *nprocs );
    void blacs_get_ ( int *ConTxt, int *what, int *val );
    void blacs_gridinit_ ( int *ConTxt, char *order, int *nprow, int *npcol );
    void blacs_gridexit_ ( int *ConTxt );
    void blacs_pcoord_ ( int *ConTxt, int *nodenum, int *prow, int *pcol );
    void blacs_barrier_ ( int*, char* );
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void dsyrk_ ( const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *lda,
                  const double *beta, double *c, const int *ldc );
    void dgemm_ ( const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a,
                  const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc );
    void dgemv_ ( const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx,
                  const double *beta, double *y, const int *incy );
    void dcopy_ ( const int *n, const double *x, const int *incx, double *y, const int *incy );
    void dlacpy_ ( const char* uplo, const int* m, const int* n,const double* a, const int* lda, double* b,const int* ldb );
    void dgemm_ ( const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda,
                  const double *b, const int *ldb, const double *beta, double *c, const int *ldc );
    void dpotrf_ ( const char* uplo, const int* n, double* a, const int* lda, int* info );
    void dpotrs_ ( const char* uplo, const int* n, const int* nrhs,const double* a, const int* lda, double* b, const int* ldb, int* info );
    void dsymm_ ( const char *side, const char *uplo, const int *m, const int *n,const double *alpha, const double *a, const int *lda, const double *b,
                  const int *ldb, const double *beta, double *c, const int *ldc );
    void pdpotrf_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdpotrs_ ( char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info );
    void pdpotri_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdsymm_( char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb,
                  int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
    void pdsymv_( char *uplo, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *x, int *ix, int *jx, int *descx, int *incx,
                  double *beta, double *y, int *iy, int *jy, int *descy, int *incy );
    void dgesd2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest );
    void dgerv2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc );

}


double d_one=1.0, d_zero=0.0, d_negone=-1.0;
int DLEN_=9, i_negone=-1, i_zero=0, i_one=1, i_two=2, i_four=4; // some many used constants
int k,n,m,l, m_plus, blocksize; //dimensions of different matrices
int lld_D, Dblocks, Ddim, t_plus, lld_S, Sblocks, Sdim, lld_B, Bblocks, Bdim;
int Drows,Dcols,Srows,Scols,Brows,Bcols;
int size, *dims, * position, ICTXT2D, ICTXT1D, ICTXT1PROC, iam;
int ntests, maxiterations,datahdf5, copyC;
char *SNPdata, *phenodata;
char *filenameT, *filenameX, *TestSet,*filenameZ, *filenamepheno;;
double lambda, epsilon;

int main( int argc, char **argv ) {
    int info, i, j, pcol;
    double * Amat, * B, *Bcopy, *Dcopy, *T, *X, *Z, *D, *y, *RHS;
    int Adim, Tdim, ell_B, RHS_dim;
    int *DESCD;

    timing secs;
    double generateRhsTime    = 0.0;
    double fillSymTime        = 0.0;
    double identityTime       = 0.0;
    double augmentedTime      = 0.0;
    double solutionTime       = 0.0;

    info = MPI_Init ( &argc, &argv );
    if ( info != 0 ) {
        printf ( "Error in MPI initialisation: %d",info );
        return info;
    }

    DESCD= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCD==NULL ) {
        printf ( "unable to allocate memory for descriptor for C\n" );
        return -1;
    }

    position= ( int* ) calloc ( 2,sizeof ( int ) );
    if ( position==NULL ) {
        printf ( "unable to allocate memory for processor position coordinate\n" );
        return EXIT_FAILURE;
    }

    dims= ( int* ) calloc ( 2,sizeof ( int ) );
    if ( dims==NULL ) {
        printf ( "unable to allocate memory for grid dimensions coordinate\n" );
        return EXIT_FAILURE;
    }

    blacs_pinfo_ ( &iam,&size ); 			//determine the number of processes involved
    blacs_setup_ ( &iam,&size );
    if ( iam ==-1 ) {
        printf ( "Error in initialisation of proces grid" );
        return -1;
    }
    info=MPI_Dims_create ( size, 2, dims );			//determine the best 2D cartesian grid with the number of processes
    if ( info != 0 ) {
        printf ( "Error in MPI creation of dimensions: %d",info );
        return info;
    }
    if (*dims != *(dims+1)) {
        while (*dims * *dims > size)
            *dims -=1;
        *(dims+1)= *dims;
        printf("WARNING: %d processor(s) unused due to reformatting to a square process grid\n", size - (*dims * *dims));
    }
    blacs_get_ ( &i_negone,&i_zero,&ICTXT2D );
    /*blacs_get_ ( &i_negone,&i_zero,&ICTXT1D );
    blacs_get_ ( &i_negone,&i_zero,&ICTXT1PROC );*/

    blacs_gridinit_ ( &ICTXT2D,"R",dims, dims+1 );
    /*blacs_gridinit_ ( &ICTXT1D,"R",&size, &i_one );
    blacs_gridinit_ ( &ICTXT1PROC,"R",&i_one, &i_one );*/

    blacs_pcoord_ ( &ICTXT2D,&iam,position, position+1 );
    if ( *position ==-1 ) {
        printf ( "Error in proces grid" );
        return -1;
    }

    info=read_input ( *++argv );
    if ( info!=0 ) {
        printf ( "Something went wrong when reading input file for processor %d\n",iam );
        return -1;
    }

    blacs_barrier_ ( &ICTXT2D,"ALL" );
    if ( * ( position+1 ) ==0 && *position==0 )
        printf ( "Reading of input-file succesful\n" );

    if ( * ( position+1 ) ==0 && *position==0 ) {
        printf("\nA linear mixed model with %d observations, %d genotypes, %d random effects and %d fixed effects\n", n,k,m,l);
        printf("was analyzed using %d (%d x %d) processors\n",size,*dims,*(dims+1));
    }

    Adim=m+l;
    RHS_dim=k+l+m;

    m_plus=k+1;
    t_plus=m+1;

    Ddim=k;
    pcol= * ( position+1 );
    Dblocks= Ddim%blocksize==0 ? Ddim/blocksize : Ddim/blocksize +1;	 //define number of blocks needed to store a complete column/row of C
    Drows= ( Dblocks - *position ) % *dims == 0 ? ( Dblocks- *position ) / *dims : ( Dblocks- *position ) / *dims +1;
    Drows= Drows<1? 1 : Drows;
    Dcols= ( Dblocks - pcol ) % * ( dims+1 ) == 0 ? ( Dblocks- pcol ) / * ( dims+1 ) : ( Dblocks- pcol ) / * ( dims+1 ) +1;
    Dcols=Dcols<1? 1 : Dcols;
    lld_D=Drows*blocksize;

    descinit_ ( DESCD, &Ddim, &Ddim, &blocksize, &blocksize, &i_zero, &i_zero, &ICTXT2D, &lld_D, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix C returns info: %d\n",info );
        return info;
    }

    D = ( double* ) calloc ( Drows * blocksize * Dcols * blocksize,sizeof ( double ) );
    if ( D==NULL ) {
        printf ( "unable to allocate memory for Matrix D  (required: %ld bytes)\n", Drows * blocksize * Dcols * blocksize * sizeof ( double ) );
        return EXIT_FAILURE;
    }

    CSRdouble BT_i, B_j;

    info = set_up_BD ( DESCD, D, BT_i, B_j );

    CSRdouble Xsparse, Zsparse;

    Xsparse.loadFromFile ( filenameX );
    Zsparse.loadFromFile ( filenameZ );

    smat_t *X_smat, *Z_smat;

    X_smat = smat_new_from ( Xsparse.nrows,Xsparse.ncols,Xsparse.pRows,Xsparse.pCols,Xsparse.pData,0,0 );
    Z_smat = smat_new_from ( Zsparse.nrows,Zsparse.ncols,Zsparse.pRows,Zsparse.pCols,Zsparse.pData,0,0 );

    smat_t *Xt_smat, *Zt_smat;
    Xt_smat = smat_copy_trans ( X_smat );
    Zt_smat = smat_copy_trans ( Z_smat );

    CSRdouble Asparse;

    smat_t *XtX_smat, *XtZ_smat, *ZtZ_smat, *ZtX_smat, *lambda_smat, *ZtZlambda_smat;
    XtX_smat = smat_matmul ( Xt_smat, X_smat );
    XtZ_smat = smat_matmul ( Xt_smat, Z_smat );
    ZtX_smat = smat_copy_trans ( XtZ_smat );
    ZtZ_smat = smat_matmul ( Zt_smat,Z_smat );

    CSRdouble Imat;

    makeIdentity ( l, Imat );

    lambda_smat = smat_new_from ( Imat.nrows,Imat.ncols,Imat.pRows,Imat.pCols,Imat.pData,0,0 );

    smat_scale_diag ( lambda_smat, -lambda );

    ZtZlambda_smat = smat_add ( lambda_smat, ZtZ_smat );

    smat_to_symmetric_structure ( XtX_smat );
    smat_to_symmetric_structure ( ZtZlambda_smat );

    CSRdouble XtX_sparse, XtZ_sparse, ZtZ_sparse;

    XtX_sparse.make ( XtX_smat->m,XtX_smat->n,XtX_smat->nnz,XtX_smat->ia,XtX_smat->ja,XtX_smat->a );
    XtZ_sparse.make ( XtZ_smat->m,XtZ_smat->n,XtZ_smat->nnz,XtZ_smat->ia,XtZ_smat->ja,XtZ_smat->a );
    ZtZ_sparse.make ( ZtZlambda_smat->m,ZtZlambda_smat->n,ZtZlambda_smat->nnz,ZtZlambda_smat->ia,ZtZlambda_smat->ja,ZtZlambda_smat->a );

    if (iam==0) {
        cout << "***                                           [  t     t  ] *** " << endl;
        cout << "***                                           [ X X   X Z ] *** " << endl;
        cout << "***                                           [           ] *** " << endl;
        cout << "*** G e n e r a t i n g    m a t r i x    A = [           ] *** " << endl;
        cout << "***                                           [  t     t  ] *** " << endl;
        cout << "***                                           [ Z X   Z Z ] *** " << endl;
    }
    blacs_barrier_ ( &ICTXT2D,"ALL" );
    create2x2SymBlockMatrix ( XtX_sparse, XtZ_sparse, ZtZ_sparse, Asparse );

    blacs_barrier_ ( &ICTXT2D,"ALL" );

    double * AB_sol;
    int * DESCAB_sol;

    AB_sol=(double *) calloc(Adim * Dcols*blocksize,sizeof(double));

    make_Sij_parallel_denseB ( Asparse,BT_i,B_j,D, lld_D, AB_sol );
    blacs_barrier_ ( &ICTXT2D,"ALL" );

    blacs_barrier_ ( &ICTXT2D,"ALL" );

    pdpotrf_ ( "U",&k,D,&i_one,&i_one,DESCD,&info );
    if ( info != 0 ) {
        printf ( "Cholesky decomposition of D was unsuccessful, error returned: %d\n",info );
        return -1;
    }
    blacs_barrier_ ( &ICTXT2D,"ALL" );

    pdpotri_ ( "U",&k,D,&i_one,&i_one,DESCD,&info );

    if ( info != 0 ) {
        printf ( "Inverse of D was unsuccessful, error returned: %d\n",info );
        return -1;
    }
    blacs_barrier_(&ICTXT2D,"A");

    double* InvD_T_Block = ( double* ) calloc ( Dblocks * blocksize + Adim ,sizeof ( double ) );

    if (iam==0) {

        int pardiso_message_level = 1;

        int pardiso_mtype=-2;

        ParDiSO pardiso ( pardiso_mtype, pardiso_message_level );
        int number_of_processors = 1;
        char* var = getenv ( "OMP_NUM_THREADS" );
        if ( var != NULL )
            sscanf ( var, "%d", &number_of_processors );

        pardiso.iparm[2]  = 2;
        pardiso.iparm[3]  = number_of_processors;
        pardiso.iparm[8]  = 0;
        pardiso.iparm[11] = 1;
        pardiso.iparm[13]  = 0;
        pardiso.iparm[28]  = 0;

        pardiso.findInverseOfA ( Asparse );

        printf("Processor %d inverted matrix A\n",iam);
    }
    blacs_barrier_(&ICTXT2D,"A");

    //Asparse.writeToFile("inverse_Asparse.txt");

    double* SX= ( double* ) calloc ( Dcols * blocksize,sizeof ( double ) );
    int * DESCSX;

    DESCAB_sol= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCAB_sol==NULL ) {
        printf ( "unable to allocate memory for descriptor for AB_sol\n" );
        return -1;
    }
    descinit_ ( DESCAB_sol, &Adim, &Ddim, &Adim, &blocksize, &i_zero, &i_zero, &ICTXT2D, &Adim, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix C returns info: %d\n",info );
        return info;
    }

    DESCSX= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCSX==NULL ) {
        printf ( "unable to allocate memory for descriptor for AB_sol\n" );
        return -1;
    }
    descinit_ ( DESCSX, &i_one, &Ddim, &i_one,&blocksize, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix C returns info: %d\n",info );
        return info;
    }
    
    blacs_barrier_(&ICTXT2D,"A");

    for (i=1; i<=Adim; ++i) {
        pdsymm_ ("R","U",&i_one,&Ddim,&d_one,D,&i_one,&i_one,DESCD,AB_sol,&i,&i_one,DESCAB_sol,&d_zero,SX,&i_one,&i_one,DESCSX);
        pddot_(&Ddim,InvD_T_Block+i-1,AB_sol,&i,&i_one,DESCAB_sol,&Adim,SX,&i_one,&i_one,DESCSX,&i_one);
    }
    
    blacs_barrier_(&ICTXT2D,"A");

    if(*position == pcol) {
        for (i=0; i<Ddim; ++i) {
            if (pcol == (i/blocksize) % *dims) {
                int Dpos = i%blocksize + ((i/blocksize) / *dims) * blocksize ;
                *(InvD_T_Block + Adim +i) = *( D + Dpos + lld_D * Dpos);
            }
        }
        for ( i=0,j=0; i<Dblocks; ++i,++j ) {
            if ( j==*dims )
                j=0;
            if ( *position==j ) {
                dgesd2d_ ( &ICTXT2D,&blocksize,&i_one,InvD_T_Block + Adim + i * blocksize,&blocksize,&i_zero,&i_zero );
            }
            if ( *position==0 ) {
                dgerv2d_ ( &ICTXT2D,&blocksize,&i_one,InvD_T_Block + Adim + blocksize*i,&blocksize,&j,&j );
            }
        }
    }

    blacs_barrier_(&ICTXT2D,"A");
    if (iam ==0) {
	for(i=0;i<Adim;++i){
	  j=Asparse.pRows[i];
	  *(InvD_T_Block+i) += Asparse.pData[j];
	}
        printdense ( Adim+k,1,InvD_T_Block,"diag_inverse_C_parallel.txt" );
    }

    return 0;
}
