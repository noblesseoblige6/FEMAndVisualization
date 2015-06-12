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

void SparseMatrix::setCCS(const vector<vector<double> >& mat)
{
  int count;
  ptr.push_back(1);
  for(int i = 0; i < colSize; i++){
    count = 0;
    for(int j = 0; j < rowSize; j++){
      if(mat[j][i] != 0){
        a.push_back(mat[j][i]);
        idx.push_back(j);
        count++;
      }
    }
    ptr.push_back(ptr[i]+count);
  }
}

void SparseMatrix::setStorage(const std::vector<double>& s)
{
  copy(s.begin(), s.end(), back_inserter(a));    
}

void SparseMatrix::setPtr(const std::vector<int>& inPtr)
{
  copy(inPtr.begin(), inPtr.end(), back_inserter(ptr));    
}

void SparseMatrix::setIdx(const std::vector<int>& inIdx)
{
  copy(inIdx.begin(), inIdx.end(), back_inserter(idx));    
}
void SparseMatrix::setParam(const vector<double>& inStorage, const vector<int>& inPtr, const vector<int>& inIdx)
{
  setStorage(inStorage);
  setPtr(inPtr);
  setIdx(inIdx);
}

vector<vector<double> > SparseMatrix::transSparseXmat(const vector<vector<double> >& mat)
{
  vector<vector<double> > res(colSize, vector<double>(mat[0].size()));
  for(int i = 0; i < colSize; i++){
    for(int j = 0; j < colSize; j++){
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
      for(int k = ptr[j]-1; k < ptr[j+1]-1; k++){
        res[i][j] += a[k]*mat[i][idx[k]];
      }
    }
  }
  return res;
}

