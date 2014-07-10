#ifndef smat_h
#define smat_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>                                                                                                                                                                                                                          
#include <string>                                                                                                                                                                                                                            
#include <time.h>                                                                                                                                                                                                                            
  //#include <omp.h>                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             
using namespace std;                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
typedef struct{                                                                                                                                                                                                                              
        double re;                                                                                                                                                                                                                           
        double i;}                                                                                                                                                                                                                           
doublecomplex;                                                                                                                                                                                                                               
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             
typedef struct pardisodef {                                                                                                                                                                                                                  
    int iparm[64], msglvl, mtype, maxfct, mnum, phase, pt[64], idum;                                                                                                                                                                         
    int * perm;                                                                                                                                                                                                                              
    double ddum;                                                                                                                                                                                                                             
} pardisoStruct;                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             
typedef struct {                                                                                                                                                                                                                             
        int   m;                /*!< Dimension rows. */                                                                                                                                                                                      
        int   n;                /*!< Dimension columns. */                                                                                                                                                                                   
        int   neig_pos; /*!< Number of positive eigenvalues. */                                                                                                                                                                              
        int   neig_neg; /*!< Number of negative eigenvalues. */
        int   perturbed_pivot;  /*!< Number of perturbed pivots. */
        int   nnz;              /*!< Number of nonzero elements.  Equal to
                                 * <tt>ia[m+1]</tt>. */

        int      sym;           /*!< symmetric? */

        int      is_complex;    /*!< complex? */
        int      *ia;           /*!< Index to rows in a and ja.  Has n+1 entries (length
                                 * of last row). */
        int        *ja;         /*!< Column of each value in a. */
        double  *a;             /*!< Values. */
} smat_t;

typedef struct {
        unsigned     n;         /* Dimension. */
        double      *v;         /* Values. */
        unsigned     pre_alloc; /* Allocated dimension.  May be used for growing. */
        int          is_complex;/* Complex? */
        /* union { */
        /*      real            *v;             /\* values *\/ */
        /*      complex double  *c; */
        /* } val; */
} vec_t;



extern "C" void  smat_to_c_indexing (smat_t * mat);
extern "C" void  smat_to_fortran_indexing (smat_t * mat); 

extern "C" smat_t * smat_add (smat_t *A, smat_t *B);

extern "C" smat_t * smat_new_from_coord (int m,
                  int nnz, 
                  int*row_index,
                  int*col_index,
                  double  *r_vals,
                  double  *i_vals,
                  int    base_index);

extern "C" smat_t * smat_copy_trans (const smat_t  *A);
extern "C" void  smat_to_symmetric_structure (smat_t  *A);
extern "C" void  smat_free     (smat_t     *A);
extern "C" smat_t*  smat_new_from_cpx (int   m,   
                int   n,   
                int  *ia, 
                int  *ja, 
                double    *a,  
                int  is_symmetric,
                int  base_index);
extern "C" void smat_print (smat_t *A);
/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int *,   int *, int *,    int *, int *, 
                             double *, int *,   int *, int *,   int *, int *,
                             int    *, double *, double *, int *, double *);
extern "C" void pardiso_residual(int *, int *, double *, int *, int *, double *, 
                                 double *, double *, double *, double *);
extern "C" void pardiso_chkmatrix_z  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec_z     (int *, int *, double *, int *);
extern "C" void pardiso_printstats_z (int *, int *, double *, int *, int *, int *,
                           double *, int *);
                           
extern "C" void smat_fc_split (smat_t   *mat, int fine_size, smat_t  **A_ff,
            smat_t  **A_fc, smat_t  **A_cf, smat_t  **A_cc);

extern "C" void smat_vert_split (smat_t     *mat, int   top_size, smat_t    **top,
              smat_t    **bot);

extern "C" void smat_horiz_split (smat_t     *mat, int   top_size, smat_t    **top,
              smat_t    **bot);

extern "C" void  smat_to_full_structure (smat_t * src,  int values,  int  *pa);

extern "C" smat_t * smat_copy_trans (const smat_t  *A);

extern "C" void smat_print_mtx(smat_t * A, char * fileName);
extern "C" void smat_print_csr(smat_t * A, char* fileName);
pardisoStruct * initPardiso(smat_t * A);
extern "C" vec_t* direct_solve (smat_t * A, pardisoStruct * pPardiso, int nrhs, vec_t * rhs, vec_t * solution);
extern "C" void vec_release (vec_t  *v);
extern "C" double vec_euclid_norm (const vec_t  *v);
extern "C" vec_t  * smat_calc_residual (vec_t  *res,  const smat_t  *A,  vec_t  *x,  const vec_t  *b);
extern "C" vec_t * vec_new (int size);          



extern "C" {
  

void smat_scale_diag (smat_t  *A, double  fac);
smat_t * smat_matmul    	(const smat_t *A, 
			 	 const smat_t *B);

smat_t * smat_copy_trans (const smat_t  *A);


smat_t  *smat_matmul3_ABCt      (const smat_t  *A, 
                                 const smat_t  *B, 
                                 const smat_t  *C);

void smat_add_diag (smat_t  *A, double  delta);


        /**
         *  Contructs a matrix from an existing csr structure.  No copies will
         *  be made.
         *
         *  <P> <B>Warning:</B> If the base has to be converted, the original
         *  structure will be altered!
         *
         *  @param      m               Row dimension
         *  @param      n               Column dimension
         *  @param      ia              Indices to Rows
         *  @param      ja              Column indices
         *  @param      a               Values
         *  @param      is_symmetric    Is the matrix symmetric?
         *  @param      base_index      Is the matrix 0-based (C) or 1-based (Fortran)?
         *
         *  @returns    A new matrix structure. */


smat_t  *smat_new_from          (unsigned        m,
				 unsigned        n,
                                 int           *ia,
                                 int           *ja,
                                 double         *a,
                                 int             is_symmetric,
                                 int             base_index);
}
#endif

