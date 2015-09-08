#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <chrono>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "rng_func.h"	//Include Path generating functions
#include "opt_eval.h"	//Include optiom calculations functions
#include "mat_opp.h"	//Include matrix operations functions
#include "download.h"	//Include Quandl API access for Hist data

using namespace std;

double CumIRate(vector <double> A, vector <double> B, double T, double dt){
	int Si = (int) A.size();
	int Brack = (int)(Si*(1 - T));
	vector<double> Rate(A.begin(), A.begin()+Brack);
	double Ret=0;
	for(size_t i = 0 ; i < Rate.size() ; i++){
		Rate[i]=A[i]-B[i];
		Ret+=Rate[i]*dt;
	}
	return Ret;
}

vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH, int T, double dt, double STRIKE){

	//Initilise counting variables
	size_t i=0, j=0;

	double Dtime = (12-DMONTH)/12;

	//Evaluate CF at time t 
	double tau = 0., n=0, tmp=0;
	vector<double> Phi;
	while(i<X[0].size()){
			//Exotic to be priced (i.e. option/spread/barrier SEE opt_eval.cpp for possible options)
		//BUTTERFLY SPREAD
			// n = -2*opt_put(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),Dtime)
				// +1*opt_put(X[1][i][T/dt - 1],STRIKE-0.1,CumIRate(X[3][i],X[4][i], Dtime, dt),Dtime)
				// +1*opt_put(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime, dt),Dtime);
		//IN-OUT Parity
			// n = barrier_call(X[1][i],1.4,0,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1)
				// n=barrier_call(X[1][i],1.3,1,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1);
			// n = barrier_call(X[1][i],1.2,0,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1);
				// +barrier_call(X[1][i],1.1,1,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1);
		//Call
			// n = opt_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1);
		//Dig Put
			n= opt_dig_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i], X[4][i], Dtime, dt), Dtime);
		//Construct Barrier with Ret clause
			// n = opt_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime, dt),1)
			// 	- opt_call(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime, dt),1)
				// - 0.1*opt_dig_call(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime, dt),1);
		//PRDC
			// n = PRDC(X[0][i],X[3][i],X[4][i],100,Dtime);
			Phi.push_back(n);
			tau+=dt;
		j=0;
		i++;
	}

	i=0;

	//Extract Variables at time t
	vector<double> tmp_vec;
	vector<double> tmp_phi;
	vector< vector<double> > tmp_mat;
	vector< vector<double> > phi_T;

	i=0;
	while(i<X.size()){
		j=0;
		tmp_phi.clear();
		while(j<X[0].size()){
			tmp_vec.push_back(X[i][j][(int)(((T/dt)/12)*DMONTH)-1]);
		 	j++;
		}
		tmp_mat.push_back(tmp_vec);
		tmp_vec.clear();
		i++;
	}

	phi_T.push_back(Phi);

	vector< vector<double> > MATRIX = transpose(tmp_mat);
	vector< vector<double> > MATRIXMultD = MatMult (tmp_mat,MATRIX);
	vector< vector<double> > MATRIXMultN = MatMult(tmp_mat ,transpose(phi_T));
	vector< vector<double> > Beta = MatAdd(MatMult(MatInv(MATRIXMultD), MATRIXMultN),MatConstMult(GenIdent((int) MATRIXMultN.size()),1e-5));

	vector <double> BetaVec;
	for(size_t i = 0; i < (Beta).size(); i++){
		for(size_t j = 0; j < (Beta)[0].size(); j++){
			tmp = (Beta)[i][j];
		}
		BetaVec.push_back(tmp);
	}

	return BetaVec;
}

vector < vector < vector < double > > >  BetaGen(int Path, int Buck, vector <double> Buck_vec, int Poly, 
												bool OUT, vector <double> Hist, int T, double dt, double STRIKE,
												bool TF){

	string strikeS = "tmp/strike.dat";
	ofstream SK;
	SK.open(strikeS.c_str());

	SK<<STRIKE<<endl;

	//Timing Start
	clock_t start;
	double duration;
	start = clock();


	vector< vector < vector <double> > > X =  PathGen(time(NULL), Path, T, dt, Hist[0]);
	vector< vector < vector <double> > > XPand = ExpPolyReg(X[0],Poly,1);
	XPand.push_back(X[1]);
	XPand.push_back(X[2]);

	//Declare variables for Bucketing
	int min, max;
	vector<double> v, Bucket;
	
	vector < vector < vector < double > > > RetVec;
	vector < vector < double > > BUCKETVEC;

	for(double DMONTH = 1; DMONTH<=12 ; DMONTH++){

		//Extract values for asset price from training paths at time t. 
		//Add to vector v
		for(size_t i=0; i<XPand[0].size(); i++){
			v.push_back(XPand[1][i][(int)((XPand[0][0].size()/12)*DMONTH)]);
		}

		//Setting auto buckets equal number in each
		vector<double> tmp_v= v;
		sort(tmp_v.begin(),tmp_v.end());
		if(TF == 0){
			for (int i = 1; i<Buck+1 ; i++){
				Buck_vec.push_back(tmp_v[(i*(int)(v.size()/(Buck+1)))]);
			}
		}

		//Find minimum and maximum values of the asset price and add buffer
		max = (int) *max_element(v.begin(),v.end()) + 1;
		min = (*min_element(v.begin(),v.end()) > 0.02)? 0.02: *min_element(v.begin(),v.end());

		//Add min, max & user elements to Bucket
		//sort in ascending order.
		Bucket.push_back(min);
		Bucket.push_back(max);
		for(int i = 0; i<Buck ; i++){
			Bucket.push_back(Buck_vec[i]);
		}

		sort(Bucket.begin(),Bucket.end());

		//SPLIT TO BUCKETS
		//And perform individual regressions, finding beta coefficents
		//and applying to function
		//Output to BucketFit.dat
		//**************NOT OPTIMISED******************

		bool YN;
		vector <double> tmp1;
		vector < vector <double> > tmp2;
		vector < vector < vector <double> > > tmp3;
		vector <double> Btmp;
		for(size_t j = 1 ; j<Bucket.size() ; j++){
			for(size_t a = 0 ; a<XPand.size(); a++){
				for( size_t i = 0 ; i<v.size() ; i++ ){
					YN = bool(Bucket[j-1]<v[i] && v[i]<=Bucket[j]);
					if(YN==1){
						for (size_t b = 0 ; b<XPand[0][0].size(); b++){
							tmp1.push_back(XPand[a][i][b]);
						}
						tmp2.push_back(tmp1);
						tmp1.clear();
					}
				}
				tmp3.push_back(tmp2);
				tmp2.clear();
			}
			//APPLY REGRESSION TO INDIVIDUAL BUCKETS
			Btmp = Reg(tmp3, DMONTH, T, dt, STRIKE);
			Btmp.insert(Btmp.begin(), Bucket[j]);
			Btmp.insert(Btmp.begin(), Bucket[j-1]);
			BUCKETVEC.push_back(Btmp);
			Btmp.clear();
			tmp3.clear();
		}
		RetVec.push_back(BUCKETVEC);
		Bucket.clear();
		Buck_vec.clear();
		v.clear();
		tmp_v.clear();
		BUCKETVEC.clear();

	}
		
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to calculate Beta coefficents: "<<duration<<"s"<<endl;

	XPand.clear();

	return RetVec;
}

vector <vector <vector <double> > > Val(double seed, double VPath, int T, double dt, double Init){
	clock_t start;
	double duration;
	start = clock();


	vector <vector <vector <double> > > ValP = PathGen(seed, VPath, T, dt, Init);


	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate valuation paths: "<<duration<<"s"<<endl;

	return ValP;

}
