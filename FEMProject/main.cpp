#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

enum LOCAL_COORD{
	BLEFT = 0,
	BRIGHT = 1,
	TRIGHT = 2,
	TLEFT = 3
};

const int ElemNum = 2;
const int NodeNum = 6;
const int NodeNumPerElem = 4;

const LOCAL_COORD NodeOrder[4] = {BLEFT, BRIGHT, TRIGHT, TLEFT};
const int NodeOrderIndx[2][4] = {{0, 1, 5, 4},{2, 3, 6, 5}};

const double LocalCoord[4][2] = {{-1 , -1},{1 , -1},{1 , 1},{-1 , 1}};
const double GlobalCoord[6][2] = {{0, 0},{1, 0},{2, 0},{0, 1},{1, 1},{2, 1}};

const double GausseXi[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
const double GausseEta[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
const double GausseOmega[2] = {1.0, 1.0};

//@comment Variables for the elastic matrix
const int E = 1;
const double Nu = 0.3;
const double lambda = Nu*E/(1-Nu*Nu);
const double Eta = E/2*(1+Nu);
const double EMat[3][3] = {
	{lambda+2*Eta, lambda      , 0},
	{lambda      , lambda+2*Eta, 0},
	{0           , 0           , Eta}
};


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

double findDetermJ(int elemIndx)
{
	double sum = 0;
	double first, second;
	for(int i = 0; i < NodeNumPerElem; i++){
		int idx = NodeOrderIndx[elemIndx][i];
		first = partialXi(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1])*GlobalCoord[idx][0]*partialEta(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1])*GlobalCoord[idx][1];
		second = partialEta(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1])*GlobalCoord[idx][0]*partialXi(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1])*GlobalCoord[idx][1];
		sum += (first - second);
	}
	return sum;
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

void findStarainDispMatrix(vector< vector<double> >& strainDispMat, vector< vector<double> >& strainDispMatT)
{
	for(int i = 0; i < NodeNumPerElem; i++){
		strainDispMat[0][i] = partialXi(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1]); 
		strainDispMat[0][i+1] = 0;

		strainDispMat[1][i] = 0; 
		strainDispMat[1][i+1] = partialEta(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1]);

		strainDispMat[2][i] = partialEta(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1]); 
		strainDispMat[2][i+1] = partialXi(NodeOrder[i], LocalCoord[i][0], LocalCoord[i][1]);
	}
	transposeMat(strainDispMat, strainDispMatT);
}

void f(double xi, double eta, int elemIdx, vector<vector<double> >& BTEB)
{
	vector< vector<double> >  B(3, vector<double>(8));
	vector< vector<double> >  BT(8, vector<double>(3));
	vector< vector<double> > EB(3, vector<double>(8));

	findStarainDispMatrix(B, BT);
	double j = findDetermJ(elemIdx);

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
			BTEB[i][j] *= j;
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
			double oemga = GausseOmega[i]*GausseOmega[j];
			f(GausseXi[i], GausseEta[j], elemIndx, BTEB);
			//@comment As n = 2, omega is 1. So ommit the multiply operation.
			//mulMat(BTEB, omega)
			addMat(res, BTEB);
		}
	}
}

void findStiffnessMatrix(vector< vector<double> >& K)
{

}

int main()
{
	vector< vector<double> > K(16, vector<double>(16));
	findStiffnessMatrix(K);
	return 0;
}