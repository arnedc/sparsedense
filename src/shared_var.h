#ifndef Shared_Var_h
#define Shared_Var_h

#define MPI_COMM_WORLD ((MPI_Comm)0x44000000)
#define MPI_DISTRIBUTE_CYCLIC   122
#define MPI_ORDER_FORTRAN        57
#  define MPI_ORDER_C             56
#define MPI_MODE_CREATE              1  /* ADIO_CREATE */ 
#define MPI_MODE_WRONLY              4  /* ADIO_WRONLY  */
typedef int MPI_Info;
#define MPI_INFO_NULL         ((MPI_Info)0x1c000000)

#ifndef MPI_FILE_DEFINED
typedef struct ADIOI_FileD *MPI_File;
#endif

#ifndef HAVE_MPI_OFFSET
typedef long long MPI_Offset;
#endif

typedef int MPI_Comm;

typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT

typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
    
} MPI_Status;

class CSRdouble;

extern "C"{
  void blacs_barrier_ ( int*, char* );
  int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
  int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
}

void printdense ( int m, int n, double *mat, char *filename );
int set_up_BD ( int* DESCD, double* Dmat, CSRdouble& BT_i, CSRdouble& B_j, CSRdouble& Btsparse ) ;
int read_input ( char* filename ) ;
int make_Sij_sparse_parallel (CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double* T_ij, int lld_Tij );
int make_Sij_parallel_denseB(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T, double* AB_sol) ;
void dense2CSR ( double *mat, int m, int n, CSRdouble& A );
void dense2CSR_sub ( double *mat, int m, int n, int lld_mat, CSRdouble& A, int startrow, int startcol );
void CSR2dense ( CSRdouble& matrix, double* T_ij ) ;
void CSR2dense_lld ( CSRdouble& matrix,double *dense, int lld_dense ) ;

void mult_colsA_colsC ( CSRdouble& A, double *B, int lld_B, int Acolstart, int Acolstop, int Ccolstart, int Ccolstop,
                        CSRdouble& C, bool trans );


extern double d_one, d_zero, d_negone;
extern int DLEN_, i_negone, i_zero, i_one; // some many used constants
extern int k,n,m, l, blocksize; //dimensions of different matrices
extern int lld_D, Dblocks, Ddim, Drows, Dcols;
extern int size, *dims, * position, ICTXT2D, iam;
extern char *filenameT, *filenameX, *filenameZ, *filenameD, *filenameC;
extern double lambda;
extern bool printD_bool, printsparseC_bool;
extern MPI_Status status;

#endif