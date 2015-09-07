#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

vector <double> diff(vector <double> Coef){
	size_t S = Coef.size();
	vector<double> NewCoef;
	for(size_t i = 1 ; i < S ; i++){
		NewCoef.push_back(i*Coef[i]);
	}
	return NewCoef;
}

double funcEval (vector <double> Coef, double x){
	size_t S = Coef.size();
	double tmp;
	for (size_t i = 0 ; i < S ; i ++){
		tmp+=Coef[i]*pow(x,i);
	}
	return tmp;
}



int main(){

	string S1 = "1.dat";
	ofstream V1;
	V1.open(S1.c_str());

	string S2 = "2.dat";
	ofstream V2;
	V2.open(S2.c_str());

	string S3 = "3.dat";
	ofstream V3;
	V3.open(S3.c_str());

	double X[8] = {-8737.34,56082.4,-152666,228429,-202856,106900,-30947.9,3796.78};
	vector <double> V(X,X + sizeof X / sizeof X[0]);

	for(double i  = 0 ; i < 2 ; i+=0.01){
		V1 <<funcEval(V,i)<<endl;
		V2 <<funcEval(diff(V),i)<<endl;
		V3 <<abs(funcEval(diff(diff(V)),i))<<endl;
	}

	return 0;

}