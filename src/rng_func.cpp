#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "download.h"

using namespace std;

vector <vector <double> > rand_gsl(int seed, int path, int T, double dt){
	double DRIFT=0;
	const gsl_rng_type * Q = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Q);
	gsl_rng_set(r, seed);

	gsl_rng_env_setup();

	double t = 0;
	int i = 0, k = 1;
	// vector <double> A;
	// vector <vector <double> > B;
	// vector <vector <vector <double> > > C;
	vector <double> A1;
	vector <vector <double> > B1;
	B1.clear();

	double tmp;

	//output in order PATH, TIME
	while(i<path){
		A1.push_back(100);
		while(t<T-dt){
			tmp = DRIFT*dt+(gsl_rng_uniform_pos(r) - 0.5);
			A1.push_back(A1[k-1]+(log(A1[k-1])*tmp));
			// A.push_back(tmp);
			t+=dt;
			k++;
		}
		t=0;
		k=1;
		// C.push_back(B);
		// B.clear();
		B1.push_back(A1);
		A1.clear();
		i++;
	}

	gsl_rng_free(r);

	return B1;
}

vector <vector <double> > CorMat(vector <vector <double> > Data){
	size_t Row = Data.size();
	size_t Col = Data[0].size();

	vector < vector <double> > Mat;
	vector <double> tmpVec;

	vector <double> tmp1;
	vector <double> tmp2;
	double * a;
	double * b;
	for (size_t i = 0; i < Row; i++){
		for (size_t j = 0; j < Row; j++){
			for (size_t k = 0; k < Col; k++){
				tmp1.push_back(Data[i][k]);
				tmp2.push_back(Data[j][k]);
			}
			a = &tmp1[0];
			b = &tmp2[0];
			tmpVec.push_back(gsl_stats_correlation(a,1,b,1,Col));
			tmp1.clear();
			tmp2.clear();
		}
		Mat.push_back(tmpVec);
		tmpVec.clear();
		
	}

	return Mat;
}

vector< vector< vector <double> > > PathGen(double seed, int PATH, double T, double dt, double tmp1){
	
	//Stochastic Interest rate 1 parameters	
	double k = 170;		//Mean reversion parameter
	double V = 0.0;	//Average interest rate
	double o = 0.05;	//Interest rate volatility

	//Stochastic Interest rate 2 parameters
	double Sk = 170;	//Mean reversion parameter
	double SV = 0.0;	//Average interest rate
	double So = 0.05;	//Interest rate volatility

	const gsl_rng_type * Q = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Q);
	gsl_rng_set(r, seed);

	const gsl_rng_type * q = gsl_rng_default;
	gsl_rng * W = gsl_rng_alloc(q);
	gsl_rng_set(W, seed*seed);

	const gsl_rng_type * O = gsl_rng_default;
	gsl_rng * S = gsl_rng_alloc(O);
	gsl_rng_set(S, seed*seed*seed);

	gsl_rng_env_setup();

	vector< vector <double> > MultOutA;
	vector< vector <double> > MultOutR1;
	vector< vector <double> > MultOutR2;
	vector< vector < vector <double> > > MultOutPLUSR;

	vector<double> R1_vec;
	vector<double> A;
	vector<double> R2_vec;
	
	double R1, R2, tmp;

	double Vol = 0.01;

	for(int COUNT = 0 ; COUNT < PATH ; COUNT++){
		R1 = 0.015;	//Initial Interest Rate T=0
		R2 = 0.02;	//Initial Volatility T=0

		//m=1;	
		
		R2_vec.push_back(R2);		
		R1_vec.push_back(R1);

		//Initial Asset Price T=0
		// tmp =1;
		tmp = tmp1;

		A.push_back(tmp);

		for( int i = 1; i<(T/dt) ; i++ ){
			//if(i==(int)((m/20)*(T/dt))){
			R1+=k*(V-R1)*dt + R1*(gsl_ran_gaussian(W,o));
			//	m++;
			//}	
			R2 += Sk*(SV-R2)*dt + R2*(gsl_ran_gaussian(S,So));
			R2=abs(R2);
			R1=abs(R1);
			R2_vec.push_back(R2);
			R1_vec.push_back(R1);
			//tmp += tmp*(R*dt+(gsl_ran_gaussian(r,0.006)));
			tmp += tmp*((R1-R2)*dt+(gsl_ran_gaussian(r,Vol)));
			A.push_back(tmp);
		}
		MultOutA.push_back(A);
		MultOutR1.push_back(R1_vec);
		MultOutR2.push_back(R2_vec);
		R1_vec.clear();
		A.clear();
		R2_vec.clear();
	}

	MultOutPLUSR.push_back(MultOutA);
	MultOutPLUSR.push_back(MultOutR1);
	MultOutPLUSR.push_back(MultOutR2);

	return MultOutPLUSR;
}