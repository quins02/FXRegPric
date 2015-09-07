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

//DECLARE FUNCTIONS
vector < vector < vector < double > > >  BetaGen();
vector < vector < vector < double > > > Val();
vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH);
void EvRiskProf();
double parab(vector<double> Coef, double X);
double hyp(vector<double> Coef, double X);
double CoefApp(vector<double> Coef, double S, double R1, double R2);
double mean(vector <double> vec);
vector <double> X_percent(vector <double> vec, double Percentile);
double CumIRate(vector <double> A, vector <double> B, double T);

//DECLARE GLOBAL VARIABLES
double STRIKE = 1.2; //Strike with Temp value
double dt = 0.01 ;	//Time step
int Path = 10;	//Number of paths with temporary value
int T = 1;	//Lenght of path
int Poly = 2; //Order of polynomial regression
bool OUT = 0; //Output paths to file flag, with temporary value
vector <double> Hist; //Holds Historical time series
vector < vector < vector < double > > > VOut; //Holds Valutation Paths
vector < vector < double > > V; //Holds Valutation Path S
vector < vector < double > > R1; //Holds Valutation Path R1
vector < vector < double > > R2; //Holds Valutation Path R2
vector < vector < vector < double > > > XMain; //Holds Beta Values
vector < vector < vector < double > > > X; //Holds Regression Path Values

int main(){
	//Timing Start
	clock_t start;
	double duration;
	start = clock();

	//Download Historical data from Quandl, to Out.csv,
	//Import as timeseries of prices as Hist.
	//Used for initial value of spot FX, and for
	//calibrating model
	Data_Download("ECB/EURUSD","8dM9KB3tW11HYvxXGCBK");
	Hist = TS_Vec("tmp/Out.csv");

	//Output aesthic title
	printf("FX Regression Risk Profiler\n"
			"***************************\n\n"
			"Current Spot Rate for EURUSD: ");
	cout<<Hist[0]<<endl;

	//********Call main functions*************
	//BetaGen:: generates regression paths and performs regression.
	//			returns value as 3-D vector in form M[Time][Bucket][Regression Coef]
	XMain = BetaGen();
	//Val:: generates valuation paths. Returns 3-D vector in form
	//		M[Asset Price/Rate 1/Rate 2][Path][Time]
	VOut = Val();
	V=VOut[0];	//Extracts only asset price from VOut
	R1=VOut[1];	//Extracts only R1 from VOut
	R2=VOut[2];	//Extracts only R2 from VOut
	//EvRiskProf:: evaluates risk profile by applying beta coefficents to valuations paths
	//			   outputs to file
	EvRiskProf();

	//Record duration
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate risk profile: "<<duration<<"s"<<endl;

	//Wait for user to manually close program
	cout<<"Press any key to exit...";
	cin.ignore();
	cin.get();

	return 0;
}

void EvRiskProf(){
	//Declare Variables necessary for valutaion
	size_t Path_Length = V[0].size();
	size_t Num_Path = V.size();
	size_t Num_Buckets = XMain[0].size();

	//Initiate Counting varaibales
	size_t Tv = 0;
	int Interval = 12;
	double tmp, tmpR1, tmpR2;

	//Out Valuation paths (Optional)
	if(OUT == 1){
		string CVANA = "tmp/VPath.dat";
		ofstream CVAA;
		CVAA.open(CVANA.c_str());

		for(size_t i = 0 ; i < V[0].size() ; i++){
			for(size_t j = 0 ; j < V.size() ; j++){
				CVAA<<V[j][i]<<"\t";
			}
			CVAA<<endl;
		}
	}


	vector < double > Val_Path;
	vector < vector < double > > Val_T;

	int i = (int) (Path_Length/Interval) ;

	for( Tv =  0 ; Tv < XMain.size()   ; Tv++ ){

		string BUCK = "tmp/BucketFit" + static_cast<ostringstream*>( &(ostringstream() << Tv+1) )->str() + string(".dat");
		ofstream BF;
		BF.open(BUCK.c_str());

	
		for( size_t j = 0 ; j < Num_Path ; j++ ){
			tmp = V[j][i];
			tmpR1 = R1[j][i];
			tmpR2 = R2[j][i];
			for ( size_t k = 0 ; k < Num_Buckets ; k++ ){
				vector < double > tmpBeta((XMain[Tv][k].begin() + 2), XMain[Tv][k].end());
				if( XMain[Tv][k][0] < tmp && tmp < XMain[Tv][k][1] ){
					Val_Path.push_back(CoefApp(tmpBeta,tmp,tmpR1,tmpR2));
					BF<<setprecision(4)<<Tv+1<<"\t"<<tmp<<"\t"<<CoefApp(tmpBeta,tmp,tmpR1,tmpR2)<<endl;
					k = Num_Buckets;
					tmpBeta.clear();
				}
				else{}
			}
		}
		Val_T.push_back(Val_Path);
		// cout<<Val_T.size()<<"\t"<<Val_T[Val_T.size()-1].size()<<endl;

		Val_Path.clear();
		BF.close();
		i += (int) (Path_Length/Interval);
	}

	string CVAN = "tmp/CVA.dat";
	ofstream CVA;
	CVA.open(CVAN.c_str());

	string CVAM = "tmp/CVAMean.dat";
	ofstream Mean;
	Mean.open(CVAM.c_str());

	string CVAPM = "tmp/CVAPMean.dat";
	ofstream PMean;
	PMean.open(CVAPM.c_str());

	for( size_t i = 0 ; i < Val_T.size() ; i++ ){
		for( size_t j = 0 ; j < Val_T[0].size() ; j++ ){
			CVA<<i+1<<"\t"<<Val_T[i][j]<<endl;
		}
		Mean<<i+1<<"\t"<<mean(Val_T[i])<<endl;
		PMean<<i+1<<"\t"<<mean(X_percent(Val_T[i],0.95))<<endl;
	}

}

vector < vector < vector < double > > >  BetaGen(){
	//Set regression options from user
	cout<<"\n----------------------------\nPricing for FX Exotics\n----------------------------\n\n";
	cout<<"Strike price: ";
	cin>>STRIKE;
	// cout<<"Month to discount to: ";
	// cin>> DMONTH;
	cout<<"Number of Paths: ";
	cin>> Path;
	cout<<"Number of Buckets: ";
	int Buck;
	cin>> Buck;
	vector<double> Buck_vec;
	cout<<"Override bucket postions (equal number)? (Y=1, N=0): ";
	bool TF;
	cin>>TF;
	if(TF == 1){
		int TEMP;
		cout<<"Bucket positons: \n";
		for(int i = 0; i<Buck ; i++){
			cin>>TEMP;
			Buck_vec.push_back(TEMP);
		}
	}
	cout<<"Order of polynomial: ";
	cin>>Poly;
	cout<<"Output all paths to file (slow)(Y=1, N=0): ";
	cin>>OUT;


	string strikeS = "tmp/strike.dat";
	ofstream SK;
	SK.open(strikeS.c_str());

	SK<<STRIKE<<endl;

	//Timing Start
	clock_t start;
	double duration;
	start = clock();


	vector< vector < vector <double> > > X =  PathGen(time(NULL), Path, T, dt, Hist[0]);
	vector< vector < vector <double> > > T = ExpPolyReg(X[0],Poly,1);
	T.push_back(X[1]);
	T.push_back(X[2]);

	//Declare variables for Bucketing
	int min, max;
	vector<double> v, Bucket;
	
	vector < vector < vector < double > > > RetVec;
	vector < vector < double > > BUCKETVEC;

	for(double DMONTH = 1; DMONTH<=12 ; DMONTH++){

		//Extract values for asset price from training paths at time t. 
		//Add to vector v
		for(size_t i=0; i<T[0].size(); i++){
			v.push_back(T[1][i][(int)((T[0][0].size()/12)*DMONTH)]);
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
		// string BUCK = "tmp/BucketFit" + static_cast<ostringstream*>( &(ostringstream() << DMONTH) )->str() + string(".dat");
		// ofstream BF;
		// BF.open(BUCK.c_str());

		//cout<<"Calculating bucketed regression coefficents...\n";

		bool YN;
		vector <double> tmp1;
		vector < vector <double> > tmp2;
		vector < vector < vector <double> > > tmp3;
		vector <double> Btmp;
		for(size_t j = 1 ; j<Bucket.size() ; j++){
			for(size_t a = 0 ; a<T.size(); a++){
				for( size_t i = 0 ; i<v.size() ; i++ ){
					YN = bool(Bucket[j-1]<v[i] && v[i]<=Bucket[j]);
					if(YN==1){
						for (size_t b = 0 ; b<T[0][0].size(); b++){
							tmp1.push_back(T[a][i][b]);
						}
						tmp2.push_back(tmp1);
						tmp1.clear();
					}
				}
				tmp3.push_back(tmp2);
				tmp2.clear();
			}
			//APPLY REGRESSION TO INDIVIDUAL BUCKETS
			Btmp = Reg(tmp3, DMONTH);
			// Generate function and output to BucketFit.dat
			// for(double s = (Bucket[j-1]) ; s < ((Bucket[j])) ; s+=0.001){
			// 	BF<<setprecision(4)<<s<<"\t"<<parab(Btmp,s)<<endl;
			// }
			Btmp.insert(Btmp.begin(), Bucket[j]);
			Btmp.insert(Btmp.begin(), Bucket[j-1]);
			BUCKETVEC.push_back(Btmp);
			//Output beta coefficents to screen
			// cout<<"Between S= "<<Bucket[j-1]<<" and "<<Bucket[j]<<" ("<<tmp3[0].size()<<" paths): "<<endl;
			// for(size_t F = 0; F<Btmp.size(); F++){
			// 	cout<<"Beta coefficent "<<j<<"["<<F<<"] at "<<DMONTH<<": "<<Btmp[F]<<endl;
			// }
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
		
		// string Full = "AllOut.dat";
		// ofstream AllOut;
		// AllOut.open(Full.c_str());

		// for (size_t i = 0 ; i < A[0].size(); i++){
		// 	for (size_t j = 0 ; j < A.size(); j++){
		// 		AllOut<<A[j][i]<<"\t";
		// 	}
		// 	AllOut<<endl;
		// }

	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to calculate Beta coefficents: "<<duration<<"s"<<endl;

	//Generate closed form solution for Call option
	//output to CFValue.dat
	// string name = "CFValue.dat";
	// ofstream CFS;
	// CFS.open(name.c_str());

	// for(double S  = 0.001 ; S<STRIKE*2 ; S+=0.001){
	// 	 CFS<<S<<"\t"<<call_CF(S, STRIKE, 0.01, ((12-DMONTH)/12), 0.06)<<endl;
	// }

	T.clear();

	return RetVec;
}

vector <vector <vector <double> > > Val(){
	clock_t start;
	double duration;
	start = clock();


	vector <vector <vector <double> > > ValP = PathGen((time(NULL)*time(NULL)),2000,T,dt, Hist[0]);
	//vector< vector <double> > XOut = ValP[0];


	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate valuation paths: "<<duration<<"s"<<endl;

	return ValP;

}

vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH){
	//Timing Start
	// clock_t start;
	// double duration;
	// start = clock();


	//Initilise counting variables
	size_t i=0, j=0;

	double Dtime = (12-DMONTH)/12;

	//Evaluate CF at time t 
	double tau = 0., n=0, tmp=0;
	vector<double> Phi;
	while(i<X[0].size()){
			//Exotic to be priced (i.e. option/spread/barrier SEE opt_eval.cpp for possible options)
		//BUTTERFLY SPREAD
			// n = -2*opt_put(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime),Dtime)
				// +1*opt_put(X[1][i][T/dt - 1],STRIKE-0.1,CumIRate(X[3][i],X[4][i], Dtime),Dtime)
				// +1*opt_put(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime),Dtime);
		//IN-OUT Parity
			// n = barrier_call(X[1][i],1.4,0,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1)
				// n=barrier_call(X[1][i],1.3,1,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1);
			// n = barrier_call(X[1][i],1.2,0,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1);
				// +barrier_call(X[1][i],1.1,1,1,STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1);
		//Call
			// n = opt_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1);
		//Dig Put
			n= opt_dig_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime),Dtime);
		//Construct Barrier with Ret clause
			// n = opt_call(X[1][i][T/dt - 1],STRIKE,CumIRate(X[3][i],X[4][i], Dtime),1)
			// 	- opt_call(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime),1)
				// - 0.1*opt_dig_call(X[1][i][T/dt - 1],STRIKE+0.1,CumIRate(X[3][i],X[4][i], Dtime),1);
		//PRDC
			// n = PRDC(X[0][i],X[3][i],X[4][i],100,Dtime);
			Phi.push_back(n);
			tau+=dt;
		//Opt<<n<<endl;
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

	// duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	// cout<<"Time for regression & training: "<<duration<<"s"<<endl;

	return BetaVec;
}

double parab(vector<double> Coef, double X){
	double Y = 0;
	for(size_t i =0 ; i<Coef.size()-2 ; i++){
		Y+=pow(X,i)*Coef[i];
	}
	if(Y<0){Y=0;}
	// if(Y>1){Y=1;}
	return Y;
}

double CoefApp(vector<double> Coef, double S, double R1, double R2){
	double Y = 0;
	for(int i =0 ; i<=Poly ; i++){
		Y+=pow(S,i)*Coef[i];
	}
	Y+=R1*Coef[Poly+1]+R2*Coef[Poly+2];
	if(Y<0){Y=0;}
	// if(Y>1){Y=1;}
	return Y;
}

double hyp(vector<double> Coef, double X){
	double TOT = 0;
	for(size_t i =0 ; i<Coef.size(); i++){
		TOT+=pow(X,i)*Coef[i];
	}
	return TOT;
}

double mean(vector <double> vec){
	int S = vec.size();
	double SUM = 0;
	for(int i = 0 ; i < S ; i++){
		SUM+=vec[i];
	}
	return SUM/S;
}

vector <double> X_percent(vector <double> vec, double Percentile){
	sort(vec.begin(),vec.end());
	int S = vec.size() - 1;
	int Bracket = S-(int)(S*Percentile);
	vector <double> Percent (vec.end()-Bracket, vec.end());
	return Percent;
}

double CumIRate(vector <double> A, vector <double> B, double T){
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
