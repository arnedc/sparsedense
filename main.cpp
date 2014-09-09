#include <stdio.h>
#include <stdlib.h>
//#include "/usr/include/openmpi-x86_64/mpi.h"
#include "src/shared_var.h"
#include "src/config.hpp"
#include "src/CSRdouble.hpp"
#include "src/IO.hpp"
#include "src/ParDiSO.hpp"
#include "src/RealMath.hpp"
#include "src/smat.h"

extern "C" {
    int MPI_Init(int *, char ***);
    int MPI_Finalize(void);
    int MPI_Dims_create(int, int, int *);
    int MPI_Barrier( MPI_Comm comm );
    int MPI_Type_create_darray(int, int, int, int [], int [], int [], int [], int,
                               MPI_Datatype, MPI_Datatype *);
    int MPI_Type_commit(MPI_Datatype *);
    int MPI_File_open(MPI_Comm, char *, int, MPI_Info, MPI_File *);
    int MPI_File_set_view(MPI_File, MPI_Offset, MPI_Datatype,
                          MPI_Datatype, char *, MPI_Info);
    int MPI_File_write_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
    int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
    void blacs_pinfo_ ( int *mypnum, int *nprocs );
    void blacs_pnum_ ( int *ConTxt, int *prow, int *pcol );
    void blacs_setup_ ( int *mypnum, int *nprocs );
    void blacs_get_ ( int *ConTxt, int *what, int *val );
    void blacs_gridinit_ ( int *ConTxt, char *order, int *nprow, int *npcol );
    void blacs_gridexit_ ( int *ConTxt );
    void blacs_pcoord_ ( int *ConTxt, int *nodenum, int *prow, int *pcol );
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void pdpotrf_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdpotri_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdsymm_( char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb,
                  int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
    void pddot_( int *n, double *dot, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void dgesd2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest );
    void dgerv2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc );
}

double d_one = 1.0, d_zero = 0.0, d_negone = -1.0;
int DLEN_=9, i_negone=-1, i_zero=0, i_one=1; // some many used constants
int k,n,m, l, blocksize; //dimensions of different matrices
int lld_D, Dblocks, Ddim, Drows, Dcols;
int size, *dims, * position, ICTXT2D, iam;
char *filenameT, *filenameX, *filenameZ, *filenameD, *filenameC;
double lambda;
bool printD_bool, printsparseC_bool;
MPI_Status status;

int main(int argc, char **argv) {
    int info, i, j, pcol, Adim;
    double *D;
    int *DESCD;
    CSRdouble BT_i, B_j, Xsparse, Zsparse, Btsparse;
    
    BT_i.allocate(0,0,0);
    B_j.allocate(0,0,0);
    Xsparse.allocate(0,0,0);
    Zsparse.allocate(0,0,0);
    Btsparse.allocate(0,0,0);

    //Initialise MPI and some MPI-variables
    info = MPI_Init ( &argc, &argv );
    if ( info != 0 ) {
        printf ( "Error in MPI initialisation: %d\n",info );
        return info;
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

    //BLACS is the interface used by PBLAS and ScaLAPACK on top of MPI

    blacs_pinfo_ ( &iam,&size ); 				//determine the number of processes involved
    info=MPI_Dims_create ( size, 2, dims );			//determine the best 2D cartesian grid with the number of processes
    if ( info != 0 ) {
        printf ( "Error in MPI creation of dimensions: %d\n",info );
        return info;
    }

    //Until now the code can only work with square process grids
    //So we try to get the biggest square grid possible with the number of processes involved
    if (*dims != *(dims+1)) {
        while (*dims * *dims > size)
            *dims -=1;
        *(dims+1)= *dims;
        if (iam==0)
            printf("WARNING: %d processor(s) unused due to reformatting to a square process grid\n", size - (*dims * *dims));
        size = *dims * *dims;
        //cout << "New size of process grid: " << size << endl;
    }

    blacs_get_ ( &i_negone,&i_zero,&ICTXT2D );

    //Initialisation of the BLACS process grid, which is referenced as ICTXT2D
    blacs_gridinit_ ( &ICTXT2D,"R",dims, dims+1 );

    if (iam < size) {

        //The rank (iam) of the process is mapped to a 2D grid: position= (process row, process column)
        blacs_pcoord_ ( &ICTXT2D,&iam,position, position+1 );
        if ( *position ==-1 ) {
            printf ( "Error in proces grid\n" );
            return -1;
        }

        //Filenames, dimensions of all matrices and other important variables are read in as global variables (see src/readinput.cpp)
        info=read_input ( *++argv );
        if ( info!=0 ) {
            printf ( "Something went wrong when reading input file for processor %d\n",iam );
            return -1;
        }

        //blacs_barrier is used to stop any process of going beyond this point before all processes have made it up to this point.
        blacs_barrier_ ( &ICTXT2D,"ALL" );
        if ( * ( position+1 ) ==0 && *position==0 )
            printf ( "Reading of input-file succesful\n" );

        if ( * ( position+1 ) ==0 && *position==0 ) {
            printf("\nA linear mixed model with %d observations, %d genotypes, %d random effects and %d fixed effects\n", n,k,m,l);
            printf("was analyzed using %d (%d x %d) processors\n",size,*dims,*(dims+1));
        }

        //Dimension of A (sparse matrix) is the number of fixed effects(m) + the sparse random effects (l)
        Adim=m+l;

        //Dimension of D (dense matrix) is the number of dense effects (k)
        Ddim=k;

        pcol= * ( position+1 );

        //Define number of blocks needed to store a complete column/row of D
        Dblocks= Ddim%blocksize==0 ? Ddim/blocksize : Ddim/blocksize +1;

        //Define the number of rowblocks needed by the current process to store its part of the dense matrix D
        Drows= ( Dblocks - *position ) % *dims == 0 ? ( Dblocks- *position ) / *dims : ( Dblocks- *position ) / *dims +1;
        Drows= Drows<1? 1 : Drows;

        //Define the number of columnblocks needed by the current process to store its part of the dense matrix D
        Dcols= ( Dblocks - pcol ) % * ( dims+1 ) == 0 ? ( Dblocks- pcol ) / * ( dims+1 ) : ( Dblocks- pcol ) / * ( dims+1 ) +1;
        Dcols=Dcols<1? 1 : Dcols;

        //Define the local leading dimension of D (keeping in mind that matrices are always stored column-wise)
        lld_D=Drows*blocksize;

        //Initialise the descriptor of the dense distributed matrix
        DESCD= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
        if ( DESCD==NULL ) {
            printf ( "unable to allocate memory for descriptor for C\n" );
            return -1;
        }

        //D with dimensions (Ddim,Ddim) is distributed over all processes in ICTXT2D, with the first element in process (0,0)
        //D is distributed into blocks of size (blocksize,blocksize), having a local leading dimension lld_D in this specific process
        descinit_ ( DESCD, &Ddim, &Ddim, &blocksize, &blocksize, &i_zero, &i_zero, &ICTXT2D, &lld_D, &info );
        if ( info!=0 ) {
            printf ( "Descriptor of matrix C returns info: %d\n",info );
            return info;
        }

        //Allocate the space necessary to store the part of D that is held into memory of this process.
        D = ( double* ) calloc ( Drows * blocksize * Dcols * blocksize,sizeof ( double ) );
        if ( D==NULL ) {
            printf ( "unable to allocate memory for Matrix D  (required: %ld bytes)\n", Drows * blocksize * Dcols * blocksize * sizeof ( double ) );
            return EXIT_FAILURE;
        }

        blacs_barrier_ ( &ICTXT2D,"ALL" );
        if (iam==0)
            printf ( "Start set up of B & D\n" );

        blacs_barrier_ ( &ICTXT2D,"ALL" );
        //set_up_BD is declared in readdist.cpp and constructs the parts of matrices B & D in each processor
        //which are necessary to create the distributed Schur complement of D
        info = set_up_BD ( DESCD, D, BT_i, B_j, Btsparse );

        //printdense(Drows*blocksize, Dcols * blocksize,D,"matrix_D.txt");

        blacs_barrier_ ( &ICTXT2D,"ALL" );
        if (iam==0)
            printf ( "Matrices B & D set up\n" );

        if(printD_bool) {

            int array_of_gsizes[2], array_of_distribs[2], array_of_dargs[2], array_of_psize[2] ;
            int buffersize;
            MPI_Datatype file_type;
            MPI_File fh;
            MPI_Status status;
            array_of_gsizes[0]=Dblocks * blocksize;
            array_of_gsizes[1]=Dblocks * blocksize;
            array_of_distribs[0]=MPI_DISTRIBUTE_CYCLIC;
            array_of_distribs[1]=MPI_DISTRIBUTE_CYCLIC;
            array_of_dargs[0]=blocksize;
            array_of_dargs[1]=blocksize;
            array_of_psize[0]=*dims;
            array_of_psize[1]=*(dims + 1);

            MPI_Type_create_darray(size,iam,2,array_of_gsizes, array_of_distribs,
                                   array_of_dargs, array_of_psize, MPI_ORDER_FORTRAN,
                                   MPI_DOUBLE, &file_type);
            MPI_Type_commit(&file_type);
            info = MPI_File_open(MPI_COMM_WORLD, filenameD,
                                 MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                 MPI_INFO_NULL, &fh);
            /*if ( ( Drows-1 ) % *(dims+1) == *position && ( Dcols-1 ) % *(dims) == pcol && Ddim%blocksize !=0 )
                buffersize=((Drows-1) * blocksize + Ddim % blocksize) * ((Dcols-1) * blocksize + Ddim % blocksize);
            else if ( ( Drows-1 ) % *(dims+1) == *position && Ddim%blocksize !=0 )
                buffersize=((Drows-1) * blocksize + Ddim % blocksize) * Dcols * blocksize;
            else if ( ( Dcols-1 ) % *(dims) == *position && Ddim%blocksize !=0 )
                buffersize=((Dcols-1) * blocksize + Ddim % blocksize) * Drows * blocksize;
            else*/
            buffersize= Dcols * Drows * blocksize * blocksize;

            MPI_File_set_view(fh, 0, MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);
            info =MPI_File_write_all(fh, D,buffersize, MPI_DOUBLE,
                                     &status);
            if(iam==0) {
                printf("Matrix D (dimension %d) is printed in file %s\n", Dblocks*blocksize,filenameD);
            }
        }



        //Now every matrix has to set up the sparse matrix A, consisting of X'X, X'Z, Z'X and Z'Z + lambda*I
        Xsparse.loadFromFile ( filenameX );
        Zsparse.loadFromFile ( filenameZ );

        smat_t *X_smat, *Z_smat;
	
	X_smat= (smat_t *) calloc(1,sizeof(smat_t));
	Z_smat= (smat_t *) calloc(1,sizeof(smat_t));
	
        X_smat = smat_new_from ( Xsparse.nrows,Xsparse.ncols,Xsparse.pRows,Xsparse.pCols,Xsparse.pData,0,0 );
        Z_smat = smat_new_from ( Zsparse.nrows,Zsparse.ncols,Zsparse.pRows,Zsparse.pCols,Zsparse.pData,0,0 );
	
        smat_t *Xt_smat, *Zt_smat;
	Xt_smat= (smat_t *) calloc(1,sizeof(smat_t));
	Zt_smat= (smat_t *) calloc(1,sizeof(smat_t));
        Xt_smat = smat_copy_trans ( X_smat );
        Zt_smat = smat_copy_trans ( Z_smat );

        CSRdouble Asparse;
        smat_t *XtX_smat, *XtZ_smat, *ZtZ_smat, *lambda_smat, *ZtZlambda_smat;
	
	XtX_smat= (smat_t *) calloc(1,sizeof(smat_t));
	XtZ_smat= (smat_t *) calloc(1,sizeof(smat_t));
	ZtZ_smat= (smat_t *) calloc(1,sizeof(smat_t));
	

        XtX_smat = smat_matmul ( Xt_smat, X_smat );
        XtZ_smat = smat_matmul ( Xt_smat, Z_smat );
        ZtZ_smat = smat_matmul ( Zt_smat,Z_smat );

        
        /*smat_free(Xt_smat);
        smat_free(Zt_smat);
        smat_free(X_smat);
        smat_free(Z_smat);*/

        CSRdouble Imat;

        makeIdentity ( l, Imat );
	
	lambda_smat= (smat_t *) calloc(1,sizeof(smat_t));

        lambda_smat = smat_new_from ( Imat.nrows,Imat.ncols,Imat.pRows,Imat.pCols,Imat.pData,0,0 );

        smat_scale_diag ( lambda_smat, -lambda );
	
	ZtZlambda_smat= (smat_t *) calloc(1,sizeof(smat_t));

        ZtZlambda_smat = smat_add ( lambda_smat, ZtZ_smat );

        /*smat_free(ZtZ_smat);
        smat_free(lambda_smat);*/
        

        smat_to_symmetric_structure ( XtX_smat );
        smat_to_symmetric_structure ( ZtZlambda_smat );

        CSRdouble XtX_sparse, XtZ_sparse, ZtZ_sparse;

        XtX_sparse.make2 ( XtX_smat->m,XtX_smat->n,XtX_smat->nnz,XtX_smat->ia,XtX_smat->ja,XtX_smat->a );
        XtZ_sparse.make2 ( XtZ_smat->m,XtZ_smat->n,XtZ_smat->nnz,XtZ_smat->ia,XtZ_smat->ja,XtZ_smat->a );
        ZtZ_sparse.make2 ( ZtZlambda_smat->m,ZtZlambda_smat->n,ZtZlambda_smat->nnz,ZtZlambda_smat->ia,ZtZlambda_smat->ja,ZtZlambda_smat->a );
	
	/*smat_free(XtX_smat);
	smat_free(XtZ_smat);
	smat_free(ZtZlambda_smat);*/
	Xsparse.clear();
	Zsparse.clear();
	Imat.clear();

        if (iam==0) {
            cout << "***                                           [  t     t  ] *** " << endl;
            cout << "***                                           [ X X   X Z ] *** " << endl;
            cout << "***                                           [           ] *** " << endl;
            cout << "*** G e n e r a t i n g    m a t r i x    A = [           ] *** " << endl;
            cout << "***                                           [  t     t  ] *** " << endl;
            cout << "***                                           [ Z X   Z Z ] *** " << endl;
        }

        //Sparse matrix A only contains the upper triangular part of A
        create2x2SymBlockMatrix ( XtX_sparse, XtZ_sparse, ZtZ_sparse, Asparse );
        //Asparse.writeToFile("A_sparse.csr");

        /*smat_free(XtX_smat);
        smat_free(XtZ_smat);
        smat_free(ZtZlambda_smat);
        XtX_sparse.clear();
        XtZ_sparse.clear();
        ZtZ_sparse.clear();*/

        blacs_barrier_ ( &ICTXT2D,"ALL" );

        if(printsparseC_bool) {
            CSRdouble Dmat, Dblock, Csparse;
            Dblock.nrows=Dblocks * blocksize;
            Dblock.ncols=Dblocks * blocksize;
	    Dblock.allocate(Dblocks * blocksize, Dblocks * blocksize, 0);
	    Dmat.allocate(0,0,0);
            for (i=0; i<Drows; ++i) {
                for(j=0; j<Dcols; ++j) {
		  if(iam==3){
		    printdense(2,2,D,"D3.txt");
		    printf("Start Dblock: (%d,%d)\n",i*blocksize,j * lld_D * blocksize);
		    printf("rows and cols of Dblock: %d\n",Dblock.nrows);
		    printf("Position in grid: (%d,%d)\n",*position, pcol);
		    printf("Position in D: (%d,%d)\n",( * ( dims) * i + *position ) *blocksize, ( * ( dims+1 ) * j + pcol ) *blocksize);}
		    //else{
                    dense2CSR_sub(D + i * blocksize + j * lld_D * blocksize,blocksize,blocksize,lld_D,Dblock,( * ( dims) * i + *position ) *blocksize,
                                  ( * ( dims+1 ) * j + pcol ) *blocksize);
                    if ( Dblock.nonzeros>0 ) {
                        if ( Dmat.nonzeros==0 ){
                            Dmat.make2 ( Dblock.nrows,Dblock.ncols,Dblock.nonzeros,Dblock.pRows,Dblock.pCols,Dblock.pData );
			    printf("Dmat is constructed from Dblock in process %d\n", iam);
			}
                        else {
                            Dmat.addBCSR ( Dblock );
			    printf("Dblock is added to Dmat in process %d\n", iam);
			}
                    }
		    
                    //Dblock.clear();
                }
            }
            blacs_barrier_(&ICTXT2D,"A");
            if(iam==3){
	      //Dmat.writeToFile("Dmat_3.txt");
	    }
            if ( iam!=0 ) {
                //Each process other than root sends its Dmat to the root process.
                MPI_Send ( & ( Dmat.nonzeros ),1, MPI_INT,0,iam,MPI_COMM_WORLD );
                MPI_Send ( & ( Dmat.pRows[0] ),Dmat.nrows + 1, MPI_INT,0,iam+size,MPI_COMM_WORLD );
                MPI_Send ( & ( Dmat.pCols[0] ),Dmat.nonzeros, MPI_INT,0,iam+2*size,MPI_COMM_WORLD );
                MPI_Send ( & ( Dmat.pData[0] ),Dmat.nonzeros, MPI_DOUBLE,0,iam+3*size,MPI_COMM_WORLD );
		//Dmat.clear();
		//Dblock.clear();
		printf("process %d sent its Dmat\n",iam);
		
            }
            else if (iam==0) {
                for ( i=1; i<size; ++i ) {
                    // The root process receives parts of Dmat sequentially from all processes and directly adds them together.
                    int nonzeroes, count;
                    MPI_Recv ( &nonzeroes,1,MPI_INT,i,i,MPI_COMM_WORLD,&status );
		    MPI_Get_count(&status, MPI_INT, &count);
		    printf("Process 0 received %d elements of process %d\n",count,i);
                    if(nonzeroes>0) {
			printf("Nonzeroes : %d\n ",nonzeroes);
                        Dblock.allocate ( Dblocks * blocksize,Dblocks * blocksize,nonzeroes );
                        MPI_Recv ( & ( Dblock.pRows[0] ), Dblocks * blocksize + 1, MPI_INT,i,i+size,MPI_COMM_WORLD,&status );
			MPI_Get_count(&status, MPI_INT, &count);
		    printf("Process 0 received %d elements of process %d\n",count,i);
                        MPI_Recv ( & ( Dblock.pCols[0] ),nonzeroes, MPI_INT,i,i+2*size,MPI_COMM_WORLD,&status );
			MPI_Get_count(&status, MPI_INT, &count);
		    printf("Process 0 received %d elements of process %d\n",count,i);
                        MPI_Recv ( & ( Dblock.pData[0] ),nonzeroes, MPI_DOUBLE,i,i+3*size,MPI_COMM_WORLD,&status );
			MPI_Get_count(&status, MPI_DOUBLE, &count);
		    printf("Process 0 received %d elements of process %d\n",count,i);

                        Dmat.addBCSR ( Dblock );
                    }
                    
                }
                //Dmat.writeToFile("D_sparse.csr");
                Dmat.reduceSymmetric();
                Btsparse.transposeIt(1);
                create2x2SymBlockMatrix(Asparse,Btsparse, Dmat, Csparse);
		Csparse.writeToFile(filenameC);
		//Csparse.clear();
            }
        }
        printf("process %d is OK\n",iam);
        blacs_barrier_(&ICTXT2D,"A");

        //AB_sol will contain the solution of A*X=B, distributed across the process rows. Processes in the same process row possess the same part of AB_sol
        double * AB_sol;
        int * DESCAB_sol;
        DESCAB_sol= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
        if ( DESCAB_sol==NULL ) {
            printf ( "unable to allocate memory for descriptor for AB_sol\n" );
            return -1;
        }
        //AB_sol (Adim, Ddim) is distributed across all processes in ICTXT2D starting from process (0,0) into blocks of size (Adim, blocksize)
        descinit_ ( DESCAB_sol, &Adim, &Ddim, &Adim, &blocksize, &i_zero, &i_zero, &ICTXT2D, &Adim, &info );
        if ( info!=0 ) {
            printf ( "Descriptor of matrix C returns info: %d\n",info );
            return info;
        }

        AB_sol=(double *) calloc(Adim * Dcols*blocksize,sizeof(double));

        // Each process calculates the Schur complement of the part of D at its disposal. (see src/schur.cpp)
        // The solution of A * Y = B_j is stored in AB_sol (= A^-1 * B_j)
	blacs_barrier_(&ICTXT2D,"A");
        make_Sij_parallel_denseB ( Asparse, BT_i, B_j, D, lld_D, AB_sol );
	BT_i.clear();
	B_j.clear();

        //From here on the Schur complement S of D is stored in D

        blacs_barrier_ ( &ICTXT2D,"ALL" );

        //The Schur complement is factorised (by ScaLAPACK)
        pdpotrf_ ( "U",&k,D,&i_one,&i_one,DESCD,&info );
        if ( info != 0 ) {
            printf ( "Cholesky decomposition of D was unsuccessful, error returned: %d\n",info );
            return -1;
        }

        //From here on the factorization of the Schur complement S is stored in D

        blacs_barrier_ ( &ICTXT2D,"ALL" );

        //The Schur complement is inverted (by ScaLAPACK)
        pdpotri_ ( "U",&k,D,&i_one,&i_one,DESCD,&info );
        if ( info != 0 ) {
            printf ( "Inverse of D was unsuccessful, error returned: %d\n",info );
            return -1;
        }

        //From here on the inverse of the Schur complement S is stored in D

        blacs_barrier_(&ICTXT2D,"A");

        double* InvD_T_Block = ( double* ) calloc ( Dblocks * blocksize + Adim ,sizeof ( double ) );

        //Diagonal elements of the (1,1) block of C^-1 are still distributed and here they are gathered in InvD_T_Block in the root process.
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

        //Only the root process performs a selected inversion of A.
        if (iam==0) {

            int pardiso_message_level = 1;

            int pardiso_mtype=-2;

            ParDiSO pardiso ( pardiso_mtype, pardiso_message_level );
            int number_of_processors = 1;
            char* var = getenv("OMP_NUM_THREADS");
            if(var != NULL)
                sscanf( var, "%d", &number_of_processors );
            else {
                printf("Set environment OMP_NUM_THREADS to 1");
                exit(1);
            }

            pardiso.iparm[2]  = 2;
            pardiso.iparm[3]  = number_of_processors;
            pardiso.iparm[8]  = 0;
            pardiso.iparm[11] = 1;
            pardiso.iparm[13]  = 0;
            pardiso.iparm[28]  = 0;

            //This function calculates the factorisation of A once again so this might be optimized.
            pardiso.findInverseOfA ( Asparse );

            printf("Processor %d inverted matrix A\n",iam);
        }
        blacs_barrier_(&ICTXT2D,"A");

        // To minimize memory usage, and because only the diagonal elements of the inverse are needed, Y' * S is calculated row by rowblocks
        // the diagonal element is calculates as the dot product of this row and the corresponding column of Y. (Y is solution of AY=B)
        double* YSrow= ( double* ) calloc ( Dcols * blocksize,sizeof ( double ) );
        int * DESCYSROW;
        DESCYSROW= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
        if ( DESCYSROW==NULL ) {
            printf ( "unable to allocate memory for descriptor for AB_sol\n" );
            return -1;
        }
        //YSrow (1,Ddim) is distributed across processes of ICTXT2D starting from process (0,0) into blocks of size (1,blocksize)
        descinit_ ( DESCYSROW, &i_one, &Ddim, &i_one,&blocksize, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
        if ( info!=0 ) {
            printf ( "Descriptor of matrix C returns info: %d\n",info );
            return info;
        }

        blacs_barrier_(&ICTXT2D,"A");

        //Calculating diagonal elements 1 by 1 of the (0,0)-block of C^-1.
        for (i=1; i<=Adim; ++i) {
            pdsymm_ ("R","U",&i_one,&Ddim,&d_one,D,&i_one,&i_one,DESCD,AB_sol,&i,&i_one,DESCAB_sol,&d_zero,YSrow,&i_one,&i_one,DESCYSROW);
            pddot_(&Ddim,InvD_T_Block+i-1,AB_sol,&i,&i_one,DESCAB_sol,&Adim,YSrow,&i_one,&i_one,DESCYSROW,&i_one);
            /*if(*position==1 && pcol==1)
            printf("Dot product in process (1,1) is: %g\n", *(InvD_T_Block+i-1));
            if(*position==0 && pcol==1)
            printf("Dot product in process (0,1) is: %g\n",*(InvD_T_Block+i-1));*/
        }
        blacs_barrier_(&ICTXT2D,"A");
	if(D != NULL)
	  free(D);
	D = NULL;



        //Only in the root process we add the diagonal elements of A^-1
        if (iam ==0) {
            for(i=0; i<Adim; ++i) {
                j=Asparse.pRows[i];
                *(InvD_T_Block+i) += Asparse.pData[j];
            }
            printdense ( Adim+k,1,InvD_T_Block,"diag_inverse_C_parallel.txt" );
        }
    }
    //cout << iam << " reached end before MPI_Barrier" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Finalize();

    return 0;
}
