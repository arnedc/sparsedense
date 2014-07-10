#ifndef IO_hpp
#define IO_hpp



void dumpCSR(const char* filename, int n, int* ia, int* ja, double* a);
void printCSR(int n, int nnz, int* ia, int* ja, double* a);



#endif
