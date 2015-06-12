#include"SparseMatrix.h"
#include <iostream>
#include <algorithm>
#include <iterator>
using namespace std;
SparseMatrix::SparseMatrix(){}
SparseMatrix::~SparseMatrix(){}

void SparseMatrix::setSize(int row, int col)
{
  rowSize = row;
  colSize = col;
}

void SparseMatrix::setStorage(const std::vector<double>& s)
{
  copy(s.begin(), s.end(), back_inserter(a));    
}

void SparseMatrix::setParam(const std::vector<int>& inPtr, const std::vector<int>& inIdx)
{
  copy(inPtr.begin(), inPtr.end(), back_inserter(ptr));    
  copy(inIdx.begin(), inIdx.end(), back_inserter(idx));    
}

vector<vector<double> > SparseMatrix::transSparseXmat(const vector<vector<double> >& mat)
{
  vector<vector<double> > res(colSize, vector<double>(colSize));
  for(int i = 0; i < colSize; i++){
    for(int j = 0; j < colSize; j++){
      res[i][j] = 0;
      for(int k = ptr[i]-1; k < ptr[i+1]-1; k++){
        res[i][j] += a[k]*mat[idx[k]][j];
      }
    }
  }
  return res;
}

vector<vector<double> > SparseMatrix::matXsparse(const vector<vector<double> >& mat)
{
  vector<vector<double> > res(mat.size(), vector<double>(colSize));

  for(int i = 0; i < mat.size(); i++){
    for(int j = 0; j < colSize; j++){
      res[i][j] = 0;
      for(int k = ptr[j]-1; k < ptr[j+1]-1; k++){
        res[i][j] += a[k]*mat[i][idx[k]];
      }
    }
  }
  return res;
}

