#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include "SparseMatrix.h"

using namespace std;

enum LOCAL_COORD{
  BLEFT = 0,
  BRIGHT = 1,
  TRIGHT = 2,
  TLEFT = 3
};
const int DimsNum = 2;
const int ElemNum = 2;
const int NodeNum = 6;
const int NodeNumPerElem = 4;

const LOCAL_COORD NodeOrder[4] = {BLEFT, BRIGHT, TRIGHT, TLEFT};
const int NodeOrderIndx[2][4] = {{0, 1, 5, 4},{1, 2, 5, 4}};

const double LocalCoord[4][2] = {{-1 , -1},{1 , -1},{1 , 1},{-1 , 1}};
const double GlobalCoord[6][2] = {{0, 0},{1, 0},{2, 0},{0, 1},{1, 1},{2, 1}};

const double GausseXi[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
const double GausseEta[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
const double GausseOmega[2] = {1.0, 1.0};

//@comment Variables for the elastic matrix
const int E = 1;
const double Nu = 0.3;
const double lambda = (Nu*E)/(1.0-Nu*Nu);
const double Eta = (double)E/(2.0*(1.0+Nu));
const double EMat[3][3] = {
  {lambda+2*Eta, lambda      , 0},
  {lambda      , lambda+2*Eta, 0},
  {0           , 0           , Eta}
};


void printMat(vector< vector<double> >& a)
{
  cout<<"Print Matrix"<<endl;
  for(int i = 0; i < a.size(); i++){
    for(int j = 0; j < a[i].size(); j++){
      cout<<a[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
}

void addMat(vector<vector<double> >& a, vector<vector<double> >& b)
{
  for(int i = 0; i < a.size(); i++){
    for(int j = 0; j < a[i].size(); j++){
      a[i][j] += b[i][j];
    }
  }
}

void mulMat(vector<vector<double> >& a, double c)
{
  for(int i = 0; i < a.size(); i++){
    for(int j = 0; j < a[i].size(); j++){
      a[i][j] *= c;
    }
  }
}

double shapeFunc(LOCAL_COORD lCoord, double xi, double eta)
{
  double coeffXi, coeffEta;
  switch(lCoord){
    case BLEFT:
      coeffXi = -1; coeffEta = -1;
      break;
    case BRIGHT:
      coeffXi = 1; coeffEta = -1;
      break;
    case TRIGHT:
      coeffXi = 1; coeffEta = 1;
      break;
    case TLEFT:
      coeffXi = -1; coeffEta = 1;
      break;
    default:
      break;
  }
  xi *= coeffXi; eta *= coeffEta;
  return 0.25 *(1+xi)*(1+eta);
}

double partialXi(LOCAL_COORD lCoord, double xi, double eta)
{
  double coeffXi, coeffEta;
  switch(lCoord){
    case BLEFT:
      coeffXi = -1; coeffEta = -1;
      break;
    case BRIGHT:
      coeffXi = 1; coeffEta = -1;
      break;
    case TRIGHT:
      coeffXi = 1; coeffEta = 1;
      break;
    case TLEFT:
      coeffXi = -1; coeffEta = 1;
      break;
    default:
      break;
  }
  eta *= coeffEta;
  return 0.25 *coeffXi*(1+eta);
}

double partialEta(LOCAL_COORD lCoord, double xi, double eta)
{
  double coeffXi, coeffEta;
  switch(lCoord){
    case BLEFT:
      coeffXi = -1; coeffEta = -1;
      break;
    case BRIGHT:
      coeffXi = 1; coeffEta = -1;
      break;
    case TRIGHT:
      coeffXi = 1; coeffEta = 1;
      break;
    case TLEFT:
      coeffXi = -1; coeffEta = 1;
      break;
    default:
      break;
  }
  xi *= coeffXi;
  return 0.25 *(1+xi)*coeffEta;
}

void findInvJ(int elemIndx, double xi, double eta, double& detJ, vector< vector<double> >& invJ)
{
  vector< vector<double> >  J(2, vector<double>(2));
  double pXi, pEta;
  for(int i = 0; i < NodeNumPerElem; i++){
    int idx = NodeOrderIndx[elemIndx][i];
    pXi = partialXi(NodeOrder[i], xi, eta);
    pEta = partialEta(NodeOrder[i], xi, eta);

    J[0][0] += pXi*GlobalCoord[idx][0]; 
    J[0][1] += pXi*GlobalCoord[idx][1]; 
    J[1][0] += pEta*GlobalCoord[idx][0]; 
    J[1][1] += pEta*GlobalCoord[idx][1];
  }
  detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

  invJ[0][0] = J[1][1]; invJ[0][1] = -J[0][1]; 
  invJ[1][0] = -J[1][0]; invJ[1][1] = J[0][0]; 
  mulMat(invJ, 1.0/detJ);
}

void transposeMat(vector< vector<double> >& o, vector< vector<double> >& t)
{
  //@comment Make transposed matrix
  for(int i = 0; i < o.size(); i++){
    for(int j = 0; j < o[i].size(); j++){
      t[j][i] = o[i][j];
    }
  }
}

void findStarainDispMatrix(vector< vector<double> >& strainDispMat, vector< vector<double> >& strainDispMatT, vector< vector<double> >& invJ, double xi, double eta)
{
  vector<double> localDispVec(2);
  vector<double> globalDispVec(2);
  for(int i = 0; i < NodeNumPerElem; i++){
    int idx = i*2;
    localDispVec[0] = partialXi(NodeOrder[i], xi, eta);
    localDispVec[1] = partialEta(NodeOrder[i], xi, eta);

    globalDispVec[0] = invJ[0][0]*localDispVec[0]+invJ[0][1]*localDispVec[1];
    globalDispVec[1] = invJ[1][0]*localDispVec[0]+invJ[1][1]*localDispVec[1];

    strainDispMat[0][idx] = globalDispVec[0]; 
    strainDispMat[0][idx+1] = 0;

    strainDispMat[1][idx] = 0; 
    strainDispMat[1][idx+1] = globalDispVec[1];

    strainDispMat[2][idx] = globalDispVec[1]; 
    strainDispMat[2][idx+1] = globalDispVec[0];
  }
  transposeMat(strainDispMat, strainDispMatT);
  //printMat(strainDispMat);
}

void f(double xi, double eta, int elemIdx, vector<vector<double> >& BTEB)
{
  vector< vector<double> >  B(3, vector<double>(8));
  vector< vector<double> >  BT(8, vector<double>(3));
  vector< vector<double> > EB(3, vector<double>(8));

  vector< vector<double> >  invJ(2, vector<double>(2));
  double detJ;
  findInvJ(elemIdx, xi, eta, detJ, invJ);

  findStarainDispMatrix(B, BT, invJ, xi, eta);

  //@comment [E][B]
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 8; j++){
      EB[i][j] = 0;
      for(int k = 0; k < 3; k++){
        EB[i][j] += EMat[i][k]*B[k][j];
      }
    }
  }
  //@comment [B^T][E][B]
  for(int i = 0; i < 8; i++){
    for(int j = 0; j < 8; j++){
      BTEB[i][j] = 0;
      for(int k = 0; k < 3; k++){
        BTEB[i][j] += BT[i][k]*EB[k][j];
      }
    }
  }
  //@comment [B^T][E][B]|J|
  for(int i = 0; i < 8; i++){
    for(int j = 0; j < 8; j++){
      BTEB[i][j] *= detJ;
    }
  }
}
void localStiffness(int elemIndx, vector<vector<double> >& res)
{
  vector< vector<double> > BTEB(8, vector<double>(8));
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++){
      f(GausseXi[i], GausseEta[j], elemIndx, BTEB);
      //@comment As n = 2, omega is 1. So ommit the multiply operation.
      //double oemga = GausseOmega[i]*GausseOmega[j];
      //mulMat(BTEB, omega)
      addMat(res, BTEB);
    }
  }
}

void combineStiffMat(vector<vector<vector<double> > >& Ktmp, vector<vector<double> >& Kd)
{
  int offset = 0;
  for(int i = 0; i < Ktmp.size(); i++){
    offset = i*DimsNum*NodeNumPerElem;
    for(int j = 0; j < Ktmp[i].size(); j++){
      for(int k = 0; k < Ktmp[i][j].size(); k++){
        Kd[j+offset][k+offset] = Ktmp[i][j][k];
      }
    }
  }
}

void boundaryCondition(vector<vector<double> >& K, vector<int>& boundConds)
{
  for(int i = 0; i < boundConds.size(); i++){
    for(int j = 0; j < K.size(); j++){
      if(boundConds[i] == j){K[boundConds[i]][j] = 1;}
      else{K[boundConds[i]][j] = K[j][boundConds[i]] = 0;}
    }
  }
}

void findStiffnessMatrix(vector< vector<double> >& K)
{
  vector< vector<double> > Kd(16, vector<double>(16));
  vector< vector<double> > KdA(16, vector<double>(12));
  vector<vector<vector<double> > > Ktmp(ElemNum, vector<vector<double> >(DimsNum*NodeNumPerElem, vector<double>(DimsNum*NodeNumPerElem)));
  //@comment find the stiffness matrix of each element 
  for(int i = 0; i < ElemNum; i++){
    localStiffness(i, Ktmp[i]);
  }
  combineStiffMat(Ktmp, Kd);

  // printMat(Kd);
  vector< vector<double> > A(ElemNum*DimsNum*NodeNumPerElem, vector<double>(12));
  A[0][0] = 1; A[1][1] = 1; A[2][2] = 1; A[3][3] = 1; 
  A[4][8] = 1; A[5][9] = 1; A[6][7] = 1; A[7][8] = 1; 

  A[8][2] = 1;  A[9][3] = 1;  A[10][4] = 1; A[11][5] = 1; 
  A[12][10] = 1; A[13][11] = 1; A[14][8] = 1; A[15][9] = 1; 

  SparseMatrix sMat;
  sMat.setSize(16, 12);
  sMat.setCCS(A);
  K = sMat.transSparseXmat(sMat.matXsparse(Kd));
}

int main()
{
  vector< vector<double> > K(12, vector<double>(12));
  findStiffnessMatrix(K);

  vector<int> boundCond(3);
  boundCond[0] = 0; boundCond[1] = 1; boundCond[2] = 6;
  boundaryCondition(K, boundCond);

  // for(int i = 0; i < 12; i++){for(int j = 0; j < 12; j++){cout<<setprecision(2)<<K[i][j]<<" & ";}cout<<"\\\\"<<endl;}
  for(int i = 0; i < 12; i++){for(int j = 0; j < 12; j++){cout<<setw(6)<<setprecision(2)<<K[i][j]<<" ";}cout<<endl;}
  return 0;
}
