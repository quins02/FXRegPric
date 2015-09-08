#ifndef REG_FUNC_H
#define REG_FUNC_H

#include <vector>
#include <string>

using namespace std;

vector < vector < vector < double > > >  BetaGen(int Path, int Buck, vector <double> Buck_vec, int Poly, 
												bool OUT, vector <double> Hist, int T, double dt, double STRIKE,
												bool TF);
vector < vector < vector < double > > > Val(double seed, double VPath, int T, double dt, double Init);
vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH, int T, double dt, double STRIKE);
double CumIRate(vector <double> A, vector <double> B, double T, double dt);

#endif