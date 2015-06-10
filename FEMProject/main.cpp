#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


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

double findDetermJ(int elemIndx, double xi, double eta)
{
	double sum = 0;
	vector< vector<double> >  J(2, vector<double>(2));
	J[0][0] = 0; J[0][1] = 0; J[1][0] = 0; J[1][1] = 0;
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
	sum = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	return J[0][0]*J[1][1] - J[0][1]*J[1][0];
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

void findStarainDispMatrix(vector< vector<double> >& strainDispMat, vector< vector<double> >& strainDispMatT, double xi, double eta)
{
	for(int i = 0; i < NodeNumPerElem; i++){
		int idx = i*2;
		strainDispMat[0][idx] = partialXi(NodeOrder[i], xi, eta); 
		strainDispMat[0][idx+1] = 0;

		strainDispMat[1][idx] = 0; 
		strainDispMat[1][idx+1] = partialEta(NodeOrder[i], xi, eta);

		strainDispMat[2][idx] = partialEta(NodeOrder[i], xi, eta); 
		strainDispMat[2][idx+1] = partialXi(NodeOrder[i], xi, eta);
	}
	transposeMat(strainDispMat, strainDispMatT);
	
	//printMat(strainDispMat);
}

void f(double xi, double eta, int elemIdx, vector<vector<double> >& BTEB)
{
	vector< vector<double> >  B(3, vector<double>(8));
	vector< vector<double> >  BT(8, vector<double>(3));
	vector< vector<double> > EB(3, vector<double>(8));

	findStarainDispMatrix(B, BT, xi, eta);
	double detJ = findDetermJ(elemIdx, xi, eta);
	
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

void combineStiffMat(vector<vector<vector<double>>>& Ktmp, vector<vector<double>>& Kd)
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


void findStiffnessMatrix(vector< vector<double> >& K)
{
	vector< vector<double> > Kd(16, vector<double>(16));
	vector< vector<double> > KdA(16, vector<double>(12));
	vector<vector<vector<double>>> Ktmp(ElemNum, vector<vector<double>>(DimsNum*NodeNumPerElem, vector<double>(DimsNum*NodeNumPerElem)));
	//@comment Initialize
	for(int i = 0; i < 2; i++){for(int j = 0; j < 8; j++){for(int k = 0; k < 8; k++){Ktmp[i][j][k] = 0;}}}
	//@comment find the stiffness matrix of each element 
	for(int i = 0; i < ElemNum; i++){
		localStiffness(i, Ktmp[i]);
	}
	combineStiffMat(Ktmp, Kd);

	printMat(Kd);

	vector< vector<double> > A(ElemNum*DimsNum*NodeNumPerElem, vector<double>(12));
	vector< vector<double> > AT(12, vector<double>(ElemNum*DimsNum*NodeNumPerElem));
	A[0][0] = 1; A[0][1] = 0; A[0][2] = 0; A[0][3] = 0; A[0][4] = 0; A[0][5] = 0; A[0][6] = 0; A[0][7] = 0; A[0][8] = 0; A[0][9] = 0; A[0][10] = 0;A[0][11] = 0;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 0; A[1][3] = 0; A[1][4] = 0; A[1][5] = 0; A[1][6] = 0; A[1][7] = 0; A[1][8] = 0; A[1][9] = 0; A[1][10] = 0;A[1][11] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 1; A[2][3] = 0; A[2][4] = 0; A[2][5] = 0; A[2][6] = 0; A[2][7] = 0; A[2][8] = 0; A[2][9] = 0; A[2][10] = 0;A[2][11] = 0;
	A[3][0] = 0; A[3][1] = 0; A[3][2] = 0; A[3][3] = 1; A[3][4] = 0; A[3][5] = 0; A[3][6] = 0; A[3][7] = 0; A[3][8] = 0; A[3][9] = 0; A[3][10] = 0;A[3][11] = 0;
	A[4][0] = 0; A[4][1] = 0; A[4][2] = 0; A[4][3] = 0; A[4][4] = 0; A[4][5] = 0; A[4][6] = 0; A[4][7] = 0; A[4][8] = 1; A[4][9] = 0; A[4][10] = 0;A[4][11] = 0;
	A[5][0] = 0; A[5][1] = 0; A[5][2] = 0; A[5][3] = 0; A[5][4] = 0; A[5][5] = 0; A[5][6] = 0; A[5][7] = 0; A[5][8] = 0; A[5][9] = 1; A[5][10] = 0;A[5][11] = 0;
	A[6][0] = 0; A[6][1] = 0; A[6][2] = 0; A[6][3] = 0; A[6][4] = 0; A[6][5] = 0; A[6][6] = 0; A[6][7] = 1; A[6][8] = 0; A[6][9] = 0; A[6][10] = 0;A[6][11] = 0;
	A[7][0] = 0; A[7][1] = 0; A[7][2] = 0; A[7][3] = 0; A[7][4] = 0; A[7][5] = 0; A[7][6] = 0; A[7][7] = 0; A[7][8] = 1; A[7][9] = 0; A[7][10] = 0;A[7][11] = 0;

	A[8][0] = 0; A[8][1] = 0; A[8][2] = 1; A[8][3] = 0; A[8][4] = 0; A[8][5] = 0; A[8][6] = 0; A[8][7] = 0; A[8][8] = 0; A[8][9] = 0; A[8][10] = 0;A[8][11] = 0;
	A[9][0] = 0; A[9][1] = 0; A[9][2] = 0; A[9][3] = 1; A[9][4] = 0; A[9][5] = 0; A[9][6] = 0; A[9][7] = 0; A[9][8] = 0; A[9][9] = 0; A[9][10] = 0;A[9][11] = 0;
	A[10][0] = 0; A[10][1] = 0; A[10][2] = 0; A[10][3] = 0; A[10][4] = 1; A[10][5] = 0; A[10][6] = 0; A[10][7] = 0; A[10][8] = 0; A[10][9] = 0; A[10][10] = 0;A[10][11] = 0;
	A[11][0] = 0; A[11][1] = 0; A[11][2] = 0; A[11][3] = 0; A[11][4] = 0; A[11][5] = 1; A[11][6] = 0; A[11][7] = 0; A[11][8] = 0; A[11][9] = 0; A[11][10] = 0;A[11][11] = 0;
	A[12][0] = 0; A[12][1] = 0; A[12][2] = 0; A[12][3] = 0; A[12][4] = 0; A[12][5] = 0; A[12][6] = 0; A[12][7] = 0; A[12][8] = 0; A[12][9] = 0; A[12][10] = 1;A[12][11] = 0;
	A[13][0] = 0; A[13][1] = 0; A[13][2] = 0; A[13][3] = 0; A[13][4] = 0; A[13][5] = 0; A[13][6] = 0; A[13][7] = 0; A[13][8] = 0; A[13][9] = 0; A[13][10] = 0;A[13][11] = 1;
	A[14][0] = 0; A[14][1] = 0; A[14][2] = 0; A[14][3] = 0; A[14][4] = 0; A[14][5] = 0; A[14][6] = 0; A[14][7] = 0; A[14][8] = 1; A[14][9] = 0; A[14][10] = 0;A[14][11] = 0;
	A[15][0] = 0; A[15][1] = 0; A[15][2] = 0; A[15][3] = 0; A[15][4] = 0; A[15][5] = 0; A[15][6] = 0; A[15][7] = 0; A[15][8] = 0; A[15][9] = 1; A[15][10] = 0;A[15][11] = 0;

	transposeMat(A, AT);

	//@comment KdA
	for(int i = 0; i < 16; i++){
		for(int j = 0; j < 12; j++){
			KdA[i][j] = 0;
			for(int k = 0; k < 16; k++){
				KdA[i][j] += Kd[i][k]*A[k][j];
			}
		}
	}
	//@comment A^TKdA
	for(int i = 0; i < 12; i++){
		for(int j = 0; j < 12; j++){
			for(int k = 0; k < 16; k++){
				K[i][j] += AT[i][k]*KdA[k][j];
			}
		}
	}
}

int main()
{
	vector< vector<double> > K(12, vector<double>(12));
	//@comment format the matrix
	for(int i = 0; i < 12; i++){for(int j = 0; j < 12; j++){K[i][j] = 0;}}
	findStiffnessMatrix(K);
	for(int i = 0; i < 12; i++){for(int j = 0; j < 12; j++){cout<<setprecision(2)<<K[i][j]<<" & ";}cout<<"\\\\"<<endl;}
	return 0;
}