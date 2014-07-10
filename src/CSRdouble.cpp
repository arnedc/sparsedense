// include files
// =============
#include "CSRdouble.hpp"
#include <cstring>
#include "searchingsorting.hpp"
#include <cassert>
// =============

CSRdouble::~CSRdouble()
{
    clear();
}


CSRdouble::CSRdouble()
{
    nrows          = 0;
    ncols          = 0;
    nonzeros       = 0;
    pRows          = 0;
    pCols          = 0;
    pData          = 0;
    name           = "UnNamed";
}


void CSRdouble::clear()
{
    delete[] pData;
    delete[] pRows;
    delete[] pCols;
}




void CSRdouble::allocate(int n, int m, int nzeros)
{
    // this is used to set the sparse structure mainly; here the pointer
    // to the values, pdata, is not necessarily ready (initialized)
    nrows            = n;
    ncols            = m;
    nonzeros         = nzeros;
    pRows            = new int[nrows + 1];
    pCols            = new int[nonzeros];
    pData            = new double[nonzeros];
}




void CSRdouble::make(int n, int m, int nzeros, int* prows,
                     int* pcols, double* pdata)
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


void CSRdouble::sortColumns()
{
    for (int i = 0; i < nrows; i++)
    {
      int index     = pRows[i];
      int entries   = pRows[i+1] - pRows[i];
      int* pcols    = pCols + index;
      double* pdata = pData + index;

      heapsort(entries, pcols, pdata);
    }
}





void CSRdouble::residual(double* r, double* x, double* b)
{
  multiply(x, r);

  for (int i = 0; i < nrows; i++)
  {
    r[i] = b[i] - r[i];
  }
}

void CSRdouble::transposeIt(int block_size)
{
  // these are standard for the transpose of a matrix
  int transpose_nrows     = ncols;
  int transpose_ncols     = nrows;
  int transpose_nonzeros  = nonzeros;

  int* transpose_prows    = new int[transpose_nrows + 1];
  int* transpose_pcols    = new int[transpose_nonzeros];
  double* transpose_pdata = new double[transpose_nonzeros*block_size];

  int* rowptr             = transpose_pcols;
  double* dataptr         = transpose_pdata;

  // now lets create the CSR structure of the transpose
  int i, j, index, from, to;
  vector<vector<int> > transpose_veccols(transpose_nrows);
  vector<vector<double> > transpose_vecdata(transpose_nrows);

  for (i = 0; i < nrows; i++)
  {
    from = pRows[i];
    to = pRows[i+1];
    for (index = from; index < to; index++)
    {
      j = pCols[index];
      transpose_veccols[j].push_back(i);
      vector<double>& v = transpose_vecdata[j];

      for (int k = 0; k < block_size; k++)
	{
	  v.push_back(pData[index*block_size + k]);
	}
    }
  }

  // we almost have our sparse structure constructed now;
  // all what is left is to copy it from the vector.
  transpose_prows[0] = 0;
  for (i = 0; i < transpose_nrows; i++)
  {
    int entries = int(transpose_veccols[i].size());
    memcpy(rowptr,  &transpose_veccols[i].front(), entries*sizeof(int));
    memcpy(dataptr, &transpose_vecdata[i].front(), entries*block_size*sizeof(double));
    rowptr  += entries;
    dataptr += entries*block_size;

    transpose_prows[i+1] = transpose_prows[i] + entries;
  }

  nrows          = transpose_nrows;
  ncols          = transpose_ncols;
  nonzeros       = transpose_nonzeros;
  pData          = transpose_pdata;
  pCols          = transpose_pcols;
  pRows          = transpose_prows;

  name           = "Transpose of an UnNamed";
  matrixType     = NORMAL;
}





void CSRdouble::multiply(double* x, double* y)
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

    default:
      cout << "matrix: \'" << name << "\' multiply, matrixType not set" << endl;
      break;
  }
}





void CSRdouble::multiplyS(double* x, double* b)
{
  memset(b, 0, nrows*sizeof(double));
  for (int i = 0; i < nrows; i++)
  {
    double x_i  = x[i];
    double a_ii = pData[pRows[i]];
    double sum  = a_ii*x_i; 

    for (int index = pRows[i]+1; index < pRows[i+1]; index++)
    {
      int j       = pCols[index];
      double a_ij = pData[index];

      sum        += a_ij*x[j];
      b[j]       += a_ij*x[i];
    }
    b[i]         += sum;
  }
}




void CSRdouble::multiplyN(double* x, double* y)
{
    for (int i = 0; i < nrows; i++)
    {
        double sum = 0.0;
        for (int index = pRows[i]; index < pRows[i+1]; index++)
        {
            int j = pCols[index];
            sum += pData[index] * x[j];
        }

        y[i] = sum;
    }
}



void CSRdouble::multiplyT(double* x, double* y)
{
    memset(y, 0, ncols*sizeof(double));
    for (int i = 0; i < nrows; i++)
    {
        for (int index = pRows[i]; index < pRows[i+1]; index++)
        {
            int j = pCols[index];
            y[j] += pData[index] * x[i];
        }
    }
}



void CSRdouble::writeToFile(const char* filename, ios::openmode mode) const
{
    cout << "\t---> Dumping matrix to file: " << filename << endl;

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
        fout.write((const char*)pData, sizeof(double)*nonzeros);
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



void CSRdouble::loadFromFile(const char* file, ios::openmode mode)
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


        /*cout << "nrows:    " << nrows    << "\n";
        cout << "ncols:    " << ncols    << "\n";
        cout << "nonzeros: " << nonzeros << "\n";*/

        pRows = new int[nrows+1];
        pCols = new int[nonzeros];
        pData = new double[nonzeros];


        fin.seekg(sizeof(int)*3);
        fin.read((char*)pRows, sizeof(int)*(nrows+1));

        fin.seekg(sizeof(int)*(nrows+1+3));
        fin.read((char*)pCols, sizeof(int)*nonzeros);

        fin.seekg(sizeof(int)*(nrows+1+3 + nonzeros));
        fin.read((char*)pData, sizeof(double)*nonzeros);

    }
    else
    {
        cout << " ascii" << std::endl;

        fin >> nrows;
        fin >> ncols;
        fin >> nonzeros;

        /*cout << "nrows:    " << nrows    << std::endl;
        cout << "ncols:    " << ncols    << std::endl;
        cout << "nonzeros: " << nonzeros << std::endl;*/

        pRows = new int[nrows+1];
        pCols = new int[nonzeros];
        pData = new double[nonzeros];

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

    int i0 = pRows[0];
    for (int i = 0; i < nrows+1; i++)
    {
        pRows[i] -= i0;
    }

    for (int i = 0; i < nonzeros; i++)
    {
        pCols[i] -= i0;
    }


}


void CSRdouble::loadFromFileCOO(const char* file)
{
    fstream fin(file, ios::in);
    cout << "opening file: " << file << " in mode: ";
    if (!fin.is_open())
    {
        cout << "couldn't open file ... " << file << "\n";
        exit(1);
    }


   cout << " ascii" << std::endl;

   fin >> nrows >> ncols >> nonzeros;

   cout << "nrows:    " << nrows    << std::endl;
   cout << "ncols:    " << ncols    << std::endl;
   cout << "nonzeros: " << nonzeros << std::endl;

   pRows = new int[nrows+1];
   pCols = new int[nonzeros];
   pData = new double[nonzeros];

   int i, j, i0;
   double aij;
   int index;
   vector<vector<int> > vvcols(nrows+1);
   vector<vector<double> > vvdata(nrows+1);
   for (index = 0; index < nonzeros; index++)
   {
       fin >> i >> j >> aij;
       vvcols[i].push_back(j);
       vvdata[i].push_back(aij);
   }

   if (vvcols[0].empty())
     i0 = 1;
   else
     i0 = 0;
   
   fin.close();

   index = 0;
   pRows[0] = 0;
   for (i = i0; i < nrows+i0; i++)
   {
     int entries = vvcols[i].size();
     heapsort(entries, &vvcols[i][0], &vvdata[i][0]);

     memcpy(&pData[index], &vvdata[i][0], entries*sizeof(double));
     memcpy(&pCols[index], &vvcols[i][0], entries*sizeof(int));

     index += entries;
     pRows[i+1-i0] = index;
   }



   for (i = 0; i < nrows+1; i++)
   {
       pRows[i] -= i0;
   }

   for (i = 0; i < nonzeros; i++)
   {
       pCols[i] -= i0;
   }
}


// This method fills the symmetric sparse structure
// so that the matrix is not any more in upper or lower
// triangular form.
void CSRdouble::fillSymmetric()
{
  int nonzeros;
  int  n        = this->nrows;
  int* prows    = this->pRows;
  int* pcols    = this->pCols;
  double* pdata = this->pData;

  vector<vector<double> > vA(n);
  vector<vector<int> >    vcols(n);

  int i;
  for (i = 0; i < n; i++)
  {
    for (int index = prows[i]; index < prows[i+1]; index++)
    {
      int j = pcols[index];

      vcols[i].push_back(j);
      double a_ij = pdata[index];
      vA[i].push_back(a_ij);

      // this is the j column in the i-th row; now we need to find the 
      // i-th column in the j-th row; If it is there we do nothing; if
      // not then we need to add it 
      if (i != j)
      {
        bool found = false;
        for (int k = prows[j]; k < prows[j+1]; k++)
        {
          int col = pcols[k];
          if (col == i)
          {
            found = true;
            break;
          }
        }

        if ( !found )
        {
          //cout << "The matrix is not Structurally Symmetric\n";
          vcols[j].push_back(i);
          vA[j].push_back(a_ij);
        }
      }
    }
  }

  int* ia = new int[n+1];
  ia[0]   = 0;
  for (i = 0; i < n; i++)
  {
    ia[i+1] = ia[i] + vcols[i].size(); 
  }

  nonzeros   = ia[n];
  int* ja    = new int[nonzeros];
  double* a  = new double[nonzeros];

  for (i = 0; i < n; i++)
  {
    int index = ia[i];
    int entries = vcols[i].size();
    for (int j = 0; j < entries; j++)
    {
      ja[index + j] = vcols[i][j];
      a[index + j]  = vA[i][j];
    }

    if (entries > 1)
      heapsort(entries, &ja[index], &a[index]);
  }

  delete[] pRows;
  delete[] pCols;
  delete[] pData;

  make(n, n, nonzeros, ia, ja, a);
  matrixType = NORMAL;
}






void CSRdouble::savedebug(const char* filename) const
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
