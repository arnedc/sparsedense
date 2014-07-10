#include <stdio.h>
#include <stdlib.h>
#include "shared_var.h"
#include "RealMath.hpp"
#include "CSRdouble.hpp"
#include <cassert>

extern "C" {
    int MPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm);
    int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
    int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void blacs_barrier_ ( int*, char* );
    int blacs_pnum_ ( int *ConTxt, int *prow, int *pcol );
    void pdgemm_ ( char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib,
                   int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
}

int set_up_BD ( int * DESCD, double * Dmat, CSRdouble& BT_i, CSRdouble& B_j ) {

    // Read-in of matrices X, Z and T from file (filename[X,Z,T])
    // X and Z are read in entrely by every process
    // T is read in strip by strip (number of rows in each process is at maximum = blocksize)
    // D is constructed directly in a distributed way
    // B is first assembled sparse in root process and afterwards the necessary parts
    // for constructing the distributed Schur complement are sent to each process

    FILE *fT;
    int ni, i,j, info;
    int *DESCT;
    double *Tblock, *temp;
    int nTblocks, nstrips, pTblocks, stripcols, lld_T, pcol, colcur,rowcur;

    MPI_Status status;

    CSRdouble Xtsparse, Ztsparse,XtT_sparse,ZtT_sparse,XtT_temp, ZtT_temp;

    Xtsparse.loadFromFile ( filenameX );
    Ztsparse.loadFromFile ( filenameZ );

    Xtsparse.transposeIt ( 1 );
    Ztsparse.transposeIt ( 1 );

    XtT_sparse.allocate ( m,k,0 );
    ZtT_sparse.allocate ( l,k,0 );



    pcol= * ( position+1 );

    // Matrix T is read in by strips of size (blocksize * *(dims+1), k)
    // Strips of T are read in row-wise and thus it is as if we store strips of T' (transpose) column-wise with dimensions (k, blocksize * *(dims+1))
    // However we must then also transpose the process grid to distribute T' correctly

    // number of strips in which we divide matrix T'
    nstrips= n % ( blocksize * * ( dims+1 ) ) ==0 ?  n / ( blocksize * * ( dims+1 ) ) : ( n / ( blocksize * * ( dims+1 ) ) ) +1;

    //the number of columns of T' included in each strip
    stripcols= blocksize * * ( dims+1 );

    //number of blocks necessary to store complete column of T'
    nTblocks= k%blocksize==0 ? k/blocksize : k/blocksize +1;

    //number of blocks necessary in this process to store complete column of T'
    pTblocks= ( nTblocks - *position ) % *dims == 0 ? ( nTblocks- *position ) / *dims : ( nTblocks- *position ) / *dims +1;
    pTblocks= pTblocks <1? 1:pTblocks;

    //local leading dimension of the strip of T' (different from process to process)
    lld_T=pTblocks*blocksize;

    // Initialisation of descriptor of strips of matrix T'
    DESCT= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCT==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }
    // strip of T (k,stripcols) is distributed across ICTXT2D starting in process (0,0) in blocks of size (blocksize,blocksize)
    // the local leading dimension in this process is lld_T
    descinit_ ( DESCT, &k, &stripcols, &blocksize, &blocksize, &i_zero, &i_zero, &ICTXT2D, &lld_T, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    // Allocation of memory for the strip of T' in all processes

    Tblock= ( double* ) calloc ( pTblocks*blocksize*blocksize, sizeof ( double ) );
    if ( Tblock==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*position,* ( position+1 ) );
        return -1;
    }

    // Initialisation of matrix D (all diagonal elements of D equal to lambda)
    temp=Dmat;
    for ( i=0,rowcur=0,colcur=0; i<Dblocks; ++i, ++colcur, ++rowcur ) {
        if ( rowcur==*dims ) {
            rowcur=0;
            temp += blocksize;
        }
        if ( colcur==* ( dims+1 ) ) {
            colcur=0;
            temp += blocksize*lld_D;
        }
        if ( *position==rowcur && * ( position+1 ) == colcur ) {
            for ( j=0; j<blocksize; ++j ) {
                * ( temp + j  * lld_D +j ) =lambda;
            }
            if ( i==Dblocks-1 && Ddim % blocksize != 0 ) {
                for ( j=blocksize-1; j>= Ddim % blocksize; --j ) {
                    * ( temp + j * lld_D + j ) =0.0;
                }
            }
        }

    }

    fT=fopen ( filenameT,"rb" );
    if ( fT==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }

    // Set up of matrix D and B per strip of T'

    for ( ni=0; ni<nstrips; ++ni ) {
        if ( ni==nstrips-1 ) {

            free ( Tblock );

            Tblock= ( double* ) calloc ( pTblocks*blocksize*blocksize, sizeof ( double ) );
            if ( Tblock==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*position,* ( position+1 ) );
                return -1;
            }
        }

        //Each process only reads in a part of the strip of T'
        //When k is not a multiple of blocksize, read-in of the last elements of the rows of T is tricky
        if ( ( nTblocks-1 ) % *dims == *position && k%blocksize !=0 ) {
            if ( ni==0 ) {
                info=fseek ( fT, ( long ) ( pcol * blocksize * ( k ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fT, ( long ) ( blocksize * ( * ( dims+1 )-1 ) * ( k ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<blocksize; ++i ) {
                info=fseek ( fT, ( long ) ( blocksize * *position * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
                for ( j=0; j < pTblocks-1; ++j ) {
                    fread ( Tblock + i*pTblocks*blocksize + j*blocksize,sizeof ( double ),blocksize,fT );
                    info=fseek ( fT, ( long ) ( ( ( *dims ) -1 ) * blocksize * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                        return -1;
                    }
                }
                fread ( Tblock + i*pTblocks*blocksize + j*blocksize,sizeof ( double ),k%blocksize,fT );
            }
            //Normal read-in of the strips of T from a binary file (each time blocksize elements are read in)
        } else {
            if ( ni==0 ) {
                info=fseek ( fT, ( long ) ( pcol * blocksize * ( k ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fT, ( long ) ( blocksize * ( * ( dims+1 )-1 ) * ( k ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<blocksize; ++i ) {
                info=fseek ( fT, ( long ) ( blocksize * *position * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
                for ( j=0; j < pTblocks-1; ++j ) {
                    fread ( Tblock + i*pTblocks*blocksize + j*blocksize,sizeof ( double ),blocksize,fT );
                    info=fseek ( fT, ( long ) ( ( * ( dims )-1 ) * blocksize * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                        return -1;
                    }
                }
                fread ( Tblock + i*pTblocks*blocksize + j*blocksize,sizeof ( double ),blocksize,fT );
                info=fseek ( fT, ( long ) ( ( k - blocksize * ( ( pTblocks-1 ) * *dims + *position +1 ) ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *position,pcol,info );
                    return -1;
                }
            }
        }

        blacs_barrier_ ( &ICTXT2D,"A" );

        // End of read-in

        // Matrix D is the sum of the multiplications of all strips of T' by their transpose
        // Up unitl now, the entire matrix is stored, not only upper/lower triangular, which is possible since D is symmetric
        // Be aware, that you akways have to allocate memory for the enitre matrix, even when only dealing with the upper/lower triangular part

        pdgemm_ ( "N","T",&k,&k,&stripcols,&d_one, Tblock,&i_one, &i_one,DESCT, Tblock,&i_one, &i_one,DESCT, &d_one, Dmat, &i_one, &i_one, DESCD ); //Z'Z
        //pdsyrk_ ( "U","N",&k,&stripcols,&d_one, Tblock,&i_one, &i_one,DESCT, &d_one, Dmat, &t_plus, &t_plus, DESCD );

        // Matrix B consists of X'T and Z'T, since each process only has some parts of T at its disposal,
        // we need to make sure that the correct colums of Z and X are multiplied by the correct columns of T.
        for ( i=0; i<pTblocks; ++i ) {
            XtT_temp.ncols=k;

            //This function multiplies the correct columns of X' with the blocks of T at the disposal of the process
            // The result is also stored immediately at the correct positions of X'T. (see src/tools.cpp)
            mult_colsA_colsC ( Xtsparse, Tblock+i*blocksize, lld_T, ( * ( dims+1 ) * ni + pcol ) *blocksize, blocksize,
                               ( *dims * i + *position ) *blocksize, blocksize, XtT_temp, 0 );
            if ( XtT_temp.nonzeros>0 ) {
                if ( XtT_sparse.nonzeros==0 )
                    XtT_sparse.make ( XtT_temp.nrows,XtT_temp.ncols,XtT_temp.nonzeros,XtT_temp.pRows,XtT_temp.pCols,XtT_temp.pData );
                else {
                    XtT_sparse.addBCSR ( XtT_temp );
                }
            }
        }
        //Same as above for calculating Z'T
        for ( i=0; i<pTblocks; ++i ) {
            ZtT_temp.ncols=k;
            mult_colsA_colsC ( Ztsparse, Tblock+i*blocksize, lld_T, ( * ( dims+1 ) * ni + pcol ) *blocksize, blocksize,
                               blocksize * ( *dims * i + *position ), blocksize, ZtT_temp, 0 );
            if ( ZtT_temp.nonzeros>0 ) {
                if ( ZtT_sparse.nonzeros==0 )
                    ZtT_sparse.make ( ZtT_temp.nrows,ZtT_temp.ncols,ZtT_temp.nonzeros,ZtT_temp.pRows,ZtT_temp.pCols,ZtT_temp.pData );
                else
                    ZtT_sparse.addBCSR ( ZtT_temp );
            }
        }
        blacs_barrier_ ( &ICTXT2D,"A" );
    }

    info=fclose ( fT );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }

    //Each process only has calculated some parts of B
    //All parts are collected by the root process (iam==0), which assembles B
    //Each process then receives BT_i and B_j corresponding to the D_ij available to the process
    if ( iam!=0 ) {
        //Each process other than root sends its X' * T and Z' * T to the root process.
        MPI_Ssend ( & ( XtT_sparse.nonzeros ),1, MPI_INT,0,iam,MPI_COMM_WORLD );
        MPI_Ssend ( & ( XtT_sparse.pRows[0] ),XtT_sparse.nrows + 1, MPI_INT,0,iam+size,MPI_COMM_WORLD );
        MPI_Ssend ( & ( XtT_sparse.pCols[0] ),XtT_sparse.nonzeros, MPI_INT,0,iam+2*size,MPI_COMM_WORLD );
        MPI_Ssend ( & ( XtT_sparse.pData[0] ),XtT_sparse.nonzeros, MPI_DOUBLE,0,iam+3*size,MPI_COMM_WORLD );
        MPI_Ssend ( & ( ZtT_sparse.nonzeros ),1, MPI_INT,0,iam,MPI_COMM_WORLD );
        MPI_Ssend ( & ( ZtT_sparse.pRows[0] ),ZtT_sparse.nrows + 1, MPI_INT,0,4*size + iam,MPI_COMM_WORLD );
        MPI_Ssend ( & ( ZtT_sparse.pCols[0] ),ZtT_sparse.nonzeros, MPI_INT,0,iam+ 5*size,MPI_COMM_WORLD );
        MPI_Ssend ( & ( ZtT_sparse.pData[0] ),ZtT_sparse.nonzeros, MPI_DOUBLE,0,iam+6*size,MPI_COMM_WORLD );

        // And eventually receives the necessary BT_i and B_j
        // Blocking sends are used, which is why the order of the receives is critical depending on the coordinates of the process
        int nonzeroes;
        if (*position >= pcol) {
            MPI_Recv ( &nonzeroes,1,MPI_INT,0,iam,MPI_COMM_WORLD,&status );
            BT_i.allocate ( blocksize*Drows,m+l,nonzeroes );
            MPI_Recv ( & ( BT_i.pRows[0] ),blocksize*Drows + 1, MPI_INT,0,iam + size,MPI_COMM_WORLD,&status );
            int count;
            MPI_Get_count(&status,MPI_INT,&count);
            BT_i.nrows=count-1;
            MPI_Recv ( & ( BT_i.pCols[0] ),nonzeroes, MPI_INT,0,iam+2*size,MPI_COMM_WORLD,&status );
            MPI_Recv ( & ( BT_i.pData[0] ),nonzeroes, MPI_DOUBLE,0,iam+3*size,MPI_COMM_WORLD,&status );

            MPI_Recv ( &nonzeroes,1, MPI_INT,0,iam+4*size,MPI_COMM_WORLD,&status );

            B_j.allocate ( blocksize*Dcols,m+l,nonzeroes );

            MPI_Recv ( & ( B_j.pRows[0] ),blocksize*Dcols + 1, MPI_INT,0,iam + 5*size,MPI_COMM_WORLD,&status );
            MPI_Get_count(&status,MPI_INT,&count);
            B_j.nrows=count-1;
            MPI_Recv ( & ( B_j.pCols[0] ),nonzeroes, MPI_INT,0,iam+6*size,MPI_COMM_WORLD,&status );
            MPI_Recv ( & ( B_j.pData[0] ),nonzeroes, MPI_DOUBLE,0,iam+7*size,MPI_COMM_WORLD,&status );

            //Actually BT_j is sent, so it still needs to be transposed
            B_j.transposeIt ( 1 );
        }
        else {
            MPI_Recv ( &nonzeroes,1, MPI_INT,0,iam+4*size,MPI_COMM_WORLD,&status );

            B_j.allocate ( blocksize*Dcols,m+l,nonzeroes );

            MPI_Recv ( & ( B_j.pRows[0] ),blocksize*Dcols + 1, MPI_INT,0,iam + 5*size,MPI_COMM_WORLD,&status );
            int count;
            MPI_Get_count(&status,MPI_INT,&count);
            B_j.nrows=count-1;

            MPI_Recv ( & ( B_j.pCols[0] ),nonzeroes, MPI_INT,0,iam+6*size,MPI_COMM_WORLD,&status );

            MPI_Recv ( & ( B_j.pData[0] ),nonzeroes, MPI_DOUBLE,0,iam+7*size,MPI_COMM_WORLD,&status );

            B_j.transposeIt ( 1 );

            MPI_Recv ( &nonzeroes,1,MPI_INT,0,iam,MPI_COMM_WORLD,&status );
            BT_i.allocate ( blocksize*Drows,m+l,nonzeroes );
            MPI_Recv ( & ( BT_i.pRows[0] ),blocksize*Drows + 1, MPI_INT,0,iam + size,MPI_COMM_WORLD,&status );
            MPI_Get_count(&status,MPI_INT,&count);
            BT_i.nrows=count-1;
            MPI_Recv ( & ( BT_i.pCols[0] ),nonzeroes, MPI_INT,0,iam+2*size,MPI_COMM_WORLD,&status );
            MPI_Recv ( & ( BT_i.pData[0] ),nonzeroes, MPI_DOUBLE,0,iam+3*size,MPI_COMM_WORLD,&status );
        }
    }
    else {
        for ( i=1; i<size; ++i ) {
            // The root process receives parts of X' * T and Z' * T sequentially from all processes and directly adds them together.
            int nonzeroes;
            MPI_Recv ( &nonzeroes,1,MPI_INT,i,i,MPI_COMM_WORLD,&status );
            if(nonzeroes>0) {
                XtT_temp.allocate ( m,k,nonzeroes );
                MPI_Recv ( & ( XtT_temp.pRows[0] ),m + 1, MPI_INT,i,i+size,MPI_COMM_WORLD,&status );
                MPI_Recv ( & ( XtT_temp.pCols[0] ),nonzeroes, MPI_INT,i,i+2*size,MPI_COMM_WORLD,&status );
                MPI_Recv ( & ( XtT_temp.pData[0] ),nonzeroes, MPI_DOUBLE,i,i+3*size,MPI_COMM_WORLD,&status );

                XtT_sparse.addBCSR ( XtT_temp );
            }

            MPI_Recv ( &nonzeroes,1, MPI_INT,i,i,MPI_COMM_WORLD,&status );

            if(nonzeroes>0) {
                ZtT_temp.allocate ( l,k,nonzeroes );

                MPI_Recv ( & ( ZtT_temp.pRows[0] ),l + 1, MPI_INT,i,4*size + i,MPI_COMM_WORLD,&status );
                MPI_Recv ( & ( ZtT_temp.pCols[0] ),nonzeroes, MPI_INT,i,i+ 5*size,MPI_COMM_WORLD,&status );
                MPI_Recv ( & ( ZtT_temp.pData[0] ),nonzeroes, MPI_DOUBLE,i,i+6*size,MPI_COMM_WORLD,&status );

                ZtT_sparse.addBCSR ( ZtT_temp );
            }
        }
        XtT_sparse.transposeIt ( 1 );
        ZtT_sparse.transposeIt ( 1 );

        // B' is created by concatening blocks X'T and Z'T
        CSRdouble Btsparse;
        create1x2BlockMatrix ( XtT_sparse, ZtT_sparse,Btsparse );

        // For each process row i BT_i is created which is also sent to processes in column i to become B_j.
        for ( int rowproc= *dims - 1; rowproc>= 0; --rowproc ) {
            BT_i.ncols=Btsparse.ncols;
            BT_i.nrows=0;
            BT_i.nonzeros=0;
            int Drows_rowproc;
            if (rowproc!=0) {
                Drows_rowproc= ( Dblocks - rowproc ) % *dims == 0 ? ( Dblocks- rowproc ) / *dims : ( Dblocks- rowproc ) / *dims +1;
                Drows_rowproc= Drows_rowproc<1? 1 : Drows_rowproc;
            }
            else
                Drows_rowproc=Drows;
            for ( i=0; i<Drows_rowproc; ++i ) {
                //Each process in row i can hold several blocks of contiguous rows of D for which we need the corresponding rows of B_T
                // Therefore we use the function extendrows to create BT_i (see src/tools.cpp)
                BT_i.extendrows ( Btsparse, ( i * *dims + rowproc ) * blocksize,blocksize );
            }
            for ( int colproc= ( rowproc==0 ? 1 : 0 ); colproc < * ( dims+1 ); ++colproc ) {
                int *curpos, rankproc;
                rankproc= blacs_pnum_ (&ICTXT2D, &rowproc,&colproc);

                MPI_Ssend ( & ( BT_i.nonzeros ),1, MPI_INT,rankproc,rankproc,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pRows[0] ),BT_i.nrows + 1, MPI_INT,rankproc,rankproc+size,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pCols[0] ),BT_i.nonzeros, MPI_INT,rankproc,rankproc+2*size,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pData[0] ),BT_i.nonzeros, MPI_DOUBLE,rankproc,rankproc+3*size,MPI_COMM_WORLD );

                //printf("BT_i's sent to processor %d\n",rankproc);

                rankproc= blacs_pnum_ (&ICTXT2D, &colproc,&rowproc);
                MPI_Ssend ( & ( BT_i.nonzeros ),1, MPI_INT,rankproc,rankproc+4*size,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pRows[0] ),BT_i.nrows + 1, MPI_INT,rankproc,rankproc+5*size,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pCols[0] ),BT_i.nonzeros, MPI_INT,rankproc,rankproc+6*size,MPI_COMM_WORLD );
                MPI_Ssend ( & ( BT_i.pData[0] ),BT_i.nonzeros, MPI_DOUBLE,rankproc,rankproc+7*size,MPI_COMM_WORLD );

                //printf("B_j's sent to processor %d\n",rankproc);
            }
        }
        B_j.make ( BT_i.nrows,BT_i.ncols,BT_i.nonzeros,BT_i.pRows,BT_i.pCols,BT_i.pData );
        B_j.transposeIt ( 1 );
    }
    free ( DESCT );
    free ( Tblock );
    return 0;
}
