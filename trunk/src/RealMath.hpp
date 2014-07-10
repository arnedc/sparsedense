#ifndef Real_Math_hpp
#define Real_Math_hpp

class CSRdouble;

void solveSystem(CSRdouble& A, double* X, double* B, int pardiso_mtype, int number_of_rhs);
void calculateSchurComplement(CSRdouble& Z, int pardiso_mtype, CSRdouble& S);
void errorReport(int number_of_rhs, CSRdouble& A, double* x, double* b);
void generateRhs(int number_of_rhs, int n, double* B);
void writeSolution(int number_of_solution_vectors, int n, double* X);



void create2x2SymBlockMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& T, // input
                             CSRdouble& C); // output

void create1x2BlockMatrix(CSRdouble& A, CSRdouble& B, // input
                             CSRdouble& C);  // output

void create2x2SymBlockMatrix_denseT(CSRdouble& A, CSRdouble& B, double* T, // input
                             CSRdouble& C); // output

void create2x2BlockMatrix_denseT_lldT(CSRdouble& A, CSRdouble& B_i_T, CSRdouble& B_j, double* T, int lld_T,// input
                             CSRdouble& C);

void createAugmentedMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& C,   // input
                                                       CSRdouble& K);  // output


void solveManyRhsUsingSchurComplement(CSRdouble& A, int nrhs, int pardiso_mtype);

void makeRhsCSRdoubleMatrix(int number_of_rhs, int n, double* B, CSRdouble& Bmat);
void makeIdentity(int n, CSRdouble& I);
void fillSymmetric(CSRdouble* pmatrix);
void norm(int n, double* x, double* y);

#endif
