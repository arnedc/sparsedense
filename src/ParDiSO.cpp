// include files
// =============
#include "ParDiSO.hpp"
#include "CSRdouble.hpp"
#include "CSRcomplex.hpp"
// =============



void ParDiSO::clear_(PardisoMemoryGroup memory_to_release)
{
  double ddum;
  int    idum;

  phase = int(memory_to_release);

  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &idum,
                &ddum,
                &idum,
                &idum,
                &idum,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);
}



void ParDiSO::clear()
{
  clear_(PARDISO_MEMORY_FOR_FACTORS);
}




double ParDiSO::memoryAllocated() const
{
  double peakMemorySymbolic = iparm[15];
  double permanentMemory = iparm[16] + iparm[17];

  double maxmem = peakMemorySymbolic > permanentMemory 
                ? peakMemorySymbolic : permanentMemory;

  return maxmem;
}





void ParDiSO::error_() const
{
  if (error != 0)
  {
    cout << "ParDiSO error: " << error << " --- ";
    switch (error)
    {
    case -1:
      cout << "input inconsistent\n";
      break;

    case -2:
      cout << "not enough memory\n";
      break;

    case -3:
      cout << "reordering problem\n";
      break;

    case -4:
      cout << "zero pivot, numerical factorization or iterative refinement problem\n";
      break;

    case -5:
      cout << "unclassified (internal) error\n";
      break;

    case -6:
      cout << "preordering failed (matrix type 11, 13 only)\n";
      break;

    case -7:
      cout << "diagonal matrix problem\n";
      break;

    case -8:
      cout << "32 bit integer overflow problem\n";
      break;

    case -10:
      cout << "No license file pardiso.lic found\n";
      break;

    case -11:
      cout << "License is expired.\n";
      break;

    case -12:
      cout << "Wrong username or hostname.\n";
      break;

    default:
      break;
    }
  }
}




ParDiSO::~ParDiSO()
{
  clear_(PARDISO_ALL_MEMORY);
}





ParDiSO::ParDiSO()
{
  perm = 0;
}





ParDiSO::ParDiSO(int pardiso_mtype, int pardiso_msglvl)
{
  // --------------------------------------------------------------------
  // ..  Setup ParDiSO control parameters und initialize the solvers     
  //     internal adress pointers. This is only necessary for the FIRST  
  //     call of the ParDiSO solver.                                     
  // --------------------------------------------------------------------
  mtype  = pardiso_mtype;
  msglvl = pardiso_msglvl;
  solver = 0;

  maxfct = 1;
  mnum   = 1;
  nrhs   = 1; 
  perm   = 0;

  PARDISOINIT_D(pt,  &mtype, &solver, &iparm[1], &dparm[1], &error);

  error_();
}








// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//                   D O U B L E    D R I V E R S
//
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




void ParDiSO::shiftIndices_(CSRdouble& A, int value)
{
  int i;
  for (i = 0; i <= A.nrows; i++)
  {
    A.pRows[i] += value;
  }

  for (i = 0; i < A.nonzeros; i++)
  {
    A.pCols[i] += value;
  }    
}







void ParDiSO::init(CSRdouble& A, int number_of_rhs)
{
  // --------------------------------------------------------------------
  // ..  Convert matrix from 0-based C-notation to Fortran 1-based       
  //     notation.                                                       
  // --------------------------------------------------------------------
  shiftIndices_(A, 1);


  // --------------------------------------------------------------------
  //  .. pardiso_chk_matrix(...)                                         
  //     Checks the consistency of the given matrix.                     
  //     Use this functionality only for debugging purposes              
  // --------------------------------------------------------------------
    
  PARDISOCHECK_D(&mtype, &A.nrows, A.pData, A.pRows, A.pCols, &error);

  if (error != 0) 
  {
    printf("\nERROR in consistency of matrix: %d", error);
    exit(1);
  }

  // --------------------------------------------------------------------
  // ..  Reordering and Symbolic Factorization.  This step also allocates
  //     all memory that is necessary for the factorization.             
  // --------------------------------------------------------------------

  double ddum;
  nrhs = number_of_rhs;
  
  phase = 11;
 
  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);
  if (error != 0) 
  {
    printf("\nERROR in symbolic factorization of matrix: %d", error);
    exit(1);
  }

  // --------------------------------------------------------------------
  // ..  Convert matrix from 1-based Fortran-notation to C 0-based       
  //     notation.                                                       
  // --------------------------------------------------------------------
  shiftIndices_(A, -1);
  error_();
}




void ParDiSO::factorize(CSRdouble& A)
{
  double ddum;

  // for factorization phase should be equal to 12
  phase = 22;

  shiftIndices_(A, 1);

  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);

  shiftIndices_(A, -1);
  
  if (error != 0) 
  {
    printf("\nERROR in factorization of matrix: %d", error);
    exit(1);
  }

  error_();
}




void ParDiSO::solve(CSRdouble& A, double* x, double* rhs)
{
  // --------------------------------------------------------------------
  // ..  Back substitution and iterative refinement.                     
  // --------------------------------------------------------------------
  phase = 33;

  shiftIndices_(A, 1);


  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                rhs,
                x,
                &error,
                &dparm[1]);


  shiftIndices_(A, -1);
  if (error != 0) 
  {
    printf("\nERROR in back-substitution of matrix: %d", error);
    exit(1);
  }
  error_();
}







bool ParDiSO::makeSchurComplement(CSRdouble& A, CSRdouble& S)
{
  double ddum;

  shiftIndices_(A, 1);
  //shiftIndices_(S, 1);
  

  // Check if this matrix is OK
  PARDISOCHECK_D(&mtype, 
                 &A.nrows, 
                 A.pData, 
                 A.pRows, 
                 A.pCols, 
                 &error);
 
  error_();
  phase     = 12;
  iparm[38] = S.nrows;


  // Perform symbolic analysis and numerical factorization
  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);


  S.nonzeros = int(iparm[39]);
  S.allocate(S.nrows, S.ncols, S.nonzeros);

  // calculate and store the Schur-complement
  PARDISOSCHUR_D(pt, 
                 &maxfct, 
                 &mnum, 
                 &mtype, 
                 S.pData, 
                 S.pRows, 
                 S.pCols);

  shiftIndices_(S, -1);
  shiftIndices_(A, -1);

  error_();

  int nonzeros = S.pRows[S.nrows];
  bool is_it_full = false;

  if (nonzeros == S.nrows*S.ncols)
    is_it_full = true;


  return is_it_full;
}





void ParDiSO::findInverseOfA(CSRdouble& A)
{
  double ddum;

  shiftIndices_(A, 1);
  

  // Check if this matrix is OK
  PARDISOCHECK_D(&mtype, 
                 &A.nrows, 
                 A.pData, 
                 A.pRows, 
                 A.pCols, 
                 &error);
 
  error_();
  phase     = 12;

  printf("Matrix checked by PARDISO\n");

  // Perform symbolic analysis and numerical factorization
  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);


  if(error !=0)
    printf("Error when factorizing matrix PARDISO: %d\n",error);
  
  error_();
  printf("Matrix factorized by PARDISO\n");
  
  
  phase     = -22;
  iparm[36] = 1;   // do not overwrite internal factor L
  
  // Perform symbolic analysis and numerical factorization
  PARDISOCALL_D(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &A.nrows,
                A.pData,
                A.pRows,
                A.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &ddum,
                &ddum,
                &error,
                &dparm[1]);

  
  if(error !=0)
    printf("Error when inverting matrix PARDISO: %d\n",error);
  
  error_();
  printf("Matrix inverted by PARDISO\n");
  shiftIndices_(A, -1);
}




// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//                 C O M P L E X    D R I V E R S
//
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



void ParDiSO::shiftIndices_(CSRcomplex& Z, int value)
{
  int i;
  for (i = 0; i <= Z.nrows; i++)
  {
    Z.pRows[i] += value;
  }

  for (i = 0; i < Z.nonzeros; i++)
  {
    Z.pCols[i] += value;
  }    
}







void ParDiSO::init(CSRcomplex& Z, int number_of_rhs)
{
  // --------------------------------------------------------------------
  // ..  Convert matrix from 0-based C-notation to Fortran 1-based       
  //     notation.                                                       
  // --------------------------------------------------------------------
  shiftIndices_(Z, 1);


  // --------------------------------------------------------------------
  //  .. pardiso_chk_matrix(...)                                         
  //     Checks the consistency of the given matrix.                     
  //     Use this functionality only for debugging purposes              
  // --------------------------------------------------------------------
    
  PARDISOCHECK_Z(&mtype, &Z.nrows, Z.pData, Z.pRows, Z.pCols, &error);

  if (error != 0) 
  {
    printf("\nERROR in consistency of matrix: %d", error);
    exit(1);
  }

  // --------------------------------------------------------------------
  // ..  Reordering and Symbolic Factorization.  This step also allocates
  //     all memory that is necessary for the factorization.             
  // --------------------------------------------------------------------

  complex<double> zdum;
  nrhs = number_of_rhs;
  
  phase = 11;
 
  PARDISOCALL_Z(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &Z.nrows,
                Z.pData,
                Z.pRows,
                Z.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &zdum,
                &zdum,
                &error,
                &dparm[1]);

  // --------------------------------------------------------------------
  // ..  Convert matrix from 1-based Fortran-notation to C 0-based       
  //     notation.                                                       
  // --------------------------------------------------------------------
  shiftIndices_(Z, -1);
  error_();
}




void ParDiSO::factorize(CSRcomplex& Z)
{
  complex<double> zdum;

  // for factorization phase should be equal to 12
  phase = 22;

  shiftIndices_(Z, 1);

  PARDISOCALL_Z(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &Z.nrows,
                Z.pData,
                Z.pRows,
                Z.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &zdum,
                &zdum,
                &error,
                &dparm[1]);

  shiftIndices_(Z, -1);

  error_();
}




void ParDiSO::solve(CSRcomplex& Z, complex<double>* x, complex<double>* rhs)
{
  // --------------------------------------------------------------------
  // ..  Back substitution and iterative refinement.                     
  // --------------------------------------------------------------------
  phase = 33;

  shiftIndices_(Z, 1);


  PARDISOCALL_Z(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &Z.nrows,
                Z.pData,
                Z.pRows,
                Z.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                rhs,
                x,
                &error,
                &dparm[1]);


  shiftIndices_(Z, -1);


  error_();
}







bool ParDiSO::makeSchurComplement(CSRcomplex& Z, CSRcomplex& S)
{
  complex<double> zdum;

  shiftIndices_(Z, 1);
  

  // Check if this matrix is OK
  PARDISOCHECK_Z(&mtype, 
                 &Z.nrows, 
                 Z.pData, 
                 Z.pRows, 
                 Z.pCols, 
                 &error);

  phase     = 12;
  iparm[38] = Z.nrows - S.nrows;


  // Perform symbolic analysis and numerical factorization
  PARDISOCALL_Z(pt,
                &maxfct,
                &mnum,
                &mtype,
                &phase,
                &Z.nrows,
                Z.pData,
                Z.pRows,
                Z.pCols,
                perm,
                &nrhs,
                &iparm[1],
                &msglvl,
                &zdum,
                &zdum,
                &error,
                &dparm[1]);


  S.nonzeros = int(iparm[39]);
  S.allocate(S.nrows, S.ncols, S.nonzeros);

  // calculate and store the Schur-complement
  PARDISOSCHUR_Z(pt, 
                 &maxfct, 
                 &mnum, 
                 &mtype, 
                 S.pData, 
                 S.pRows, 
                 S.pCols);

  shiftIndices_(S, -1);
  shiftIndices_(Z, -1);

  error_();

  int nonzeros = S.pRows[S.nrows];
  bool is_it_full = false;

  if (nonzeros == S.nrows*S.ncols)
    is_it_full = true;


  return is_it_full;
}



