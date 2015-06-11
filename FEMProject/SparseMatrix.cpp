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

void SparseMatrix::setRowParam(const std::vector<int>& ptr, const std::vector<int>& idx)
{
  copy(ptr.begin(), ptr.end(), back_inserter(rowPtr));    
  copy(idx.begin(), idx.end(), back_inserter(rowIdx));    
}

void SparseMatrix::setColParam(const std::vector<int>& ptr, const std::vector<int>& idx)
{
  copy(ptr.begin(), ptr.end(), back_inserter(colPtr));    
  copy(idx.begin(), idx.end(), back_inserter(colIdx));    
}

vector<vector<double> > SparseMatrix::sparseXmat(const vector<vector<double> >& mat)
{
  vector<vector<double> > res(rowPtr.size(), vector<double>(mat.size()));
  for(int i = 0; i < rowPtr.size(); i++){
    for(int j = 0; j < mat.size(); j++){
      res[i][j] = 0;
      for(int k = rowPtr[i]-1; k < rowPtr[i+1]-1; k++){
        res[i][j] += a[k]*mat[rowIdx[k]][j];
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
      for(int k = colPtr[j]-1; k < colPtr[j+1]-1; k++){
        // cout<<i<<" "<<j<<" "<<k<<endl;
        res[i][j] += a[k]*mat[i][colIdx[k]];
      }
    }
  }
  return res;
}

