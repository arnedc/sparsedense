// include files
// =============
#include "CSRcomplex.hpp"
#include "searchingsorting.hpp"
#include <cstring>
// =============

CSRcomplex::~CSRcomplex()
{
    clear();
}



CSRcomplex::CSRcomplex()
{
    nrows          = 0;
    ncols          = 0;
    nonzeros       = 0;
    pRows          = 0;
    pCols          = 0;
    pData          = 0;
    name           = "UnNamed";
}


void
CSRcomplex::clear()
{
    delete[] pData;
    delete[] pRows;
    delete[] pCols;
}




void
CSRcomplex::allocate(int n, int m, int nzeros)
{
    // this is used to set the sparse structure mainly; here the pointer
    // to the values, pdata, is not necessarily ready (initialized)
    nrows            = n;
    ncols            = m;
    nonzeros         = nzeros;
    pRows            = new int[nrows + 1];
    pCols            = new int[nonzeros];
    pData            = new complex<double>[nonzeros];
}




void
CSRcomplex::make(int n, int m, int nzeros, int* prows,
                 int* pcols, complex<double>* pdata)
{
    // this is used to set the sparse structure mainly; here the pointer
    // to the values, pdata, is not necessarily ready (initialized)
    nrows            = n;
    ncols            = m;
    nonzeros         = nzeros;
    pRows            = prows;
    pCols            = pcols;
    pData            = pdata;
    name             = "UnNamed";
}




void CSRcomplex::sortColumns()
{
    for (int i = 0; i < nrows; i++)
    {
      int index              = pRows[i];
      int entries            = pRows[i+1] - pRows[i];
      int* pcols             = pCols + index;
      complex<double>* pdata = pData + index;

      heapsort(entries, pcols, pdata);
    }
}





void CSRcomplex::residual(complex<double>* r, complex<double>* x, complex<double>* b)
{
  multiply(x, r);

  for (int i = 0; i < nrows; i++)
  {
    r[i] = b[i] - r[i];
  }
}



void CSRcomplex::multiply(complex<double>* x, complex<double>* y)
{
  switch (matrixType)
  {
    case NORMAL:
      multiplyN(x, y);
      break;

    case TRANSPOSE:
      multiplyT(x, y);
      break;

    case SYMMETRIC:
      multiplyS(x, y);
      break;

    case HERMITIAN:
      multiplyH(x, y);
      break;

    default:
      cout << "matrix: \'" << name << "\' multiply, matrixType not set" << endl;
      break;
  }
}



void CSRcomplex::multiplyS(complex<double>* x, complex<double>* b)
{
  // memset(b, 0, n*sizeof(complex<double>));
  for (int i = 0; i < nrows; i++)
  {
    b[i] = 0.0;
  }

  for (int i = 0; i < nrows; i++)
  {
    complex<double> x_i  = x[i];
    complex<double> a_ii = pData[pRows[i]];
    complex<double> sum  = a_ii*x_i; 

    for (int index = pRows[i]+1; index < pRows[i+1]; index++)
    {
      int j                 = pCols[index];
      complex<double>& a_ij = pData[index];

      sum                  += a_ij*x[j];
      b[j]                 += a_ij*x_i;
    }
    b[i]                   += sum;
  }
}


void CSRcomplex::multiplyH(complex<double>* x, complex<double>* b)
{
  // memset(b, 0, n*sizeof(complex<double>));
  for (int i = 0; i < nrows; i++)
  {
    b[i] = 0.0;
  }

  for (int i = 0; i < nrows; i++)
  {
    complex<double> x_i  = x[i];
    complex<double> a_ii = pData[pRows[i]];
    complex<double> sum  = a_ii*x_i; 

    for (int index = pRows[i]+1; index < pRows[i+1]; index++)
    {
      int j                 = pCols[index];
      complex<double>& a_ij = pData[index];

      sum                  += a_ij*x[j];
      b[j]                 += conj(a_ij)*x_i;
    }
    b[i]                   += sum;
  }
}



void CSRcomplex::multiplyN(complex<double>* x, complex<double>* y)
{
    for (int i = 0; i < nrows; i++)
    {
        complex<double> sum = 0.0;
        for (int index = pRows[i]; index < pRows[i+1]; index++)
        {
            int j = pCols[index];
            sum += pData[index] * x[j];
        }

        y[i] = sum;
    }
}



void CSRcomplex::multiplyT(complex<double>* x, complex<double>* y)
{
    memset(y, 0, ncols*sizeof(complex<double>));
    for (int i = 0; i < nrows; i++)
    {
        for (int index = pRows[i]; index < pRows[i+1]; index++)
        {
            int j = pCols[index];
            y[j] += pData[index] * x[i];
        }
    }
}



void CSRcomplex::writeToFile(const char* filename, ios::openmode mode) const
{
    fstream fout(filename, ios::out | mode);
    if (!fout.is_open())
    {
        cout << "could not open file " << filename << " for output\n";
        return;
    }

    if (mode == ios::binary)
    {
        fout.seekp(0);
        fout.write((char*)&nrows, sizeof(int));

        fout.seekp(sizeof(int));
        fout.write((char*)&ncols, sizeof(int));

        fout.seekp(sizeof(int)*2);
        fout.write((char*)&nonzeros, sizeof(int));

        fout.seekp(sizeof(int)*3);
        fout.write((const char*)pRows, sizeof(int)*(nrows+1));

        fout.seekp(sizeof(int)*(nrows+1+3));
        fout.write((const char*)pCols, sizeof(int)*nonzeros);

        fout.seekp(sizeof(int)*(nrows+1+3 + nonzeros));
        fout.write((const char*)pData, sizeof(complex<double>)*nonzeros);
        fout.close();
    }
    else
    {
        fout << nrows << "\n";
        fout << ncols << "\n";
        fout << nonzeros << "\n";

        int i;
        for (i = 0; i < nrows+1; i++)
        {
            fout << pRows[i] << "\n";
        }

        for (i = 0; i < nonzeros; i++)
        {
            fout << pCols[i] << "\n";
        }

        fout.setf(ios::scientific, ios::floatfield);
        fout.precision(16);

        for (i = 0; i < nonzeros; i++)
        {
            fout << pData[i] << "\n";
        }
    }

    fout.close();
}



void CSRcomplex::loadFromFile(const char* file, ios::openmode mode)
{
    fstream fin(file, ios::in | mode);

    cout << "opening file: " << file << " in mode: ";
    if (!fin.is_open())
    {
        cout << "couldn't open file ... " << file << "\n";
        exit(1);
    }


    if (mode == ios::binary)
    {
        cout << " binary" << std::endl;
        fin.seekg(0);
        fin.read((char*)&nrows, sizeof(int));

        fin.seekg(sizeof(int));
        fin.read((char*)&ncols, sizeof(int));

        fin.seekg(sizeof(int)*2);
        fin.read((char*)&nonzeros, sizeof(int));


        cout << "nrows:    " << nrows    << "\n";
        cout << "ncols:    " << ncols    << "\n";
        cout << "nonzeros: " << nonzeros << "\n";

        pRows = new int[nrows+1];
        pCols = new int[nonzeros];
        pData = new complex<double>[nonzeros];


        fin.seekg(sizeof(int)*3);
        fin.read((char*)pRows, sizeof(int)*(nrows+1));

        fin.seekg(sizeof(int)*(nrows+1+3));
        fin.read((char*)pCols, sizeof(int)*nonzeros);

        fin.seekg(sizeof(int)*(nrows+1+3 + nonzeros));
        fin.read((char*)pData, sizeof(complex<double>)*nonzeros);

    }
    else
    {
        cout << " ascii" << std::endl;

        fin >> nrows;
        fin >> ncols;
        fin >> nonzeros;

        cout << "nrows:    " << nrows    << "\n";
        cout << "ncols:    " << ncols    << "\n";
        cout << "nonzeros: " << nonzeros << "\n";

        pRows = new int[nrows+1];
        pCols = new int[nonzeros];
        pData = new complex<double>[nonzeros];

        int i;
        for (i = 0; i < nrows+1; i++)
        {
            fin >> pRows[i];
        }

        for (i = 0; i < nonzeros; i++)
        {
            fin >> pCols[i];
        }

        for (i = 0; i < nonzeros; i++)
        {
            fin >> pData[i];
        }
    }
    fin.close();


    for (int i = 0; i < nrows+1; i++)
    {
        pRows[i] -= pRows[0];
    }

    for (int i = 0; i < nonzeros; i++)
    {
        pCols[i] -= pRows[0];
    }
}





void CSRcomplex::savedebug(const char* filename) const
{
  fstream fout(filename, ios::out);

  int i, index, j;

  fout.setf(ios::scientific, ios::floatfield);
  fout.precision(16);

  for (i = 0; i < nrows; i++)
  {
      fout << "row #" << i << "\n";
      fout << "============\n";

      for (index = pRows[i]; index < pRows[i+1]; index++)
      {
          j = pCols[index];
          fout << setw(12) << j << setw(1) << ":" << setw(25) << pData[index] << "\n";
      }
      fout << "\n";
  }
}
