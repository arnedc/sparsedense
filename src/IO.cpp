#include "config.hpp"
#include "IO.hpp"


void dumpCSR(const char* filename, int n, int* ia, int* ja, double* a)
{
  fstream fout(filename, ios::out);
  fout << n << endl;
  fout << n << endl;
  fout << ia[n] << endl;
  
  for (int i = 0; i <= n; i++)
  {
    fout << ia[i] << endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << ja[i] << endl;
  }

  for (int i = 0; i < ia[n]; i++)
  {
    fout << a[i] << endl;
  }

  fout.close();
}




void printCSR(int n, int nnz, int* ia, int* ja, double* a)
{
  cout << "rows: " << setw(10) << n << endl;
  cout << "nnz : " << setw(10) << nnz << endl;

  if (nnz == n*n)
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        cout << setw(10) << a[i*n + j];
      }
      cout << endl;
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      for (int index = ia[i]; index < ia[i+1]; index++)
      {
        int j = ja[index];
        cout << setw(10) << "(" << i << ", " << j << ") " << a[index];
      }
      cout << endl;
    }
  }
}


