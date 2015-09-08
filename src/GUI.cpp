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
#include "reg_func.h"	//Include regression functions

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>

using namespace std;

//DECLARE FUNCTIONS
double EvRiskProf(double Poly, double dt);
double parab(vector<double> Coef, double X);
double hyp(vector<double> Coef, double X);
double CoefApp(vector<double> Coef, double S, double R1, double R2, double Poly);
double mean(vector <double> vec);
vector <double> X_percent(vector <double> vec, double Percentile);

//DECLARE GLOBAL VARIABLES
vector < vector < vector < double > > > VOut; //Holds Valutation Paths
vector < vector < double > > V; //Holds Valutation Path S
vector < vector < double > > R1; //Holds Valutation Path R1
vector < vector < double > > R2; //Holds Valutation Path R2
vector < vector < vector < double > > > XMain; //Holds Beta Values
vector < vector < vector < double > > > X; //Holds Regression Path Values

int main(int argc, char **argv){
	//Timing Start
	clock_t start;
	double duration;
	start = clock();

	double STRIKE = 1.2; //Strike with Temp value
	double dt = 0.01 ;	//Time step
	int Path = 10;	//Number of paths with temporary value
	int T = 1;	//Lenght of path
	int Poly = 2; //Order of polynomial regression
	bool OUT = 0; //Output paths to file flag, with temporary value
	vector <double> Hist; //Holds Historical time series

	//Download Historical data from Quandl, to Out.csv,
	//Import as timeseries of prices as Hist.
	//Used for initial value of spot FX, and for
	//calibrating model
	Data_Download("ECB/EURUSD","8dM9KB3tW11HYvxXGCBK");
	Hist = TS_Vec("tmp/Out.csv");

	int Buck= 1;
	vector<double> Buck_vec;
	bool TF=0;
	if(TF == 1){
		int TEMP=1;
		for(int i = 0; i<Buck ; i++){
			Buck_vec.push_back(TEMP);
		}
	}

	Fl_Window *window = new Fl_Window(680,360);
	Fl_Input *input = new Fl_Input(50, 10, 100, 30, "Strike");
	input->value("1,2");

	window->end();
  	window->show(argc, argv);

	//********Call main functions*************

	//BetaGen:: generates regression paths and performs regression.
	//			returns value as 3-D vector in form M[Time][Bucket][Regression Coef]
	XMain = BetaGen(Path, //Number of Regression Paths
					Buck, //Number of Buckets
					Buck_vec, //Optionally filled vector of bucket positions
					Poly, //Order of polynomial for regression
					OUT,	//boolean for file out
					Hist,	//vector of historic data
					T,	//Nominal total time, usually T=1
					dt,	//Time step
					STRIKE,	//Strike of Option
					TF);	//Boolean for auto filled buckets

	//Val:: generates valuation paths. Returns 3-D vector in form
	//		M[Asset Price/Rate 1/Rate 2][Path][Time]
	VOut = Val((time(NULL)*time(NULL)), //Seed of RNG
				2000,	//Number of Paths
				T,	//Nominal total time, usuallt T=1
				dt, //Time Step
				Hist[0]);	//Initial point for stoch proccess
	V=VOut[0];	//Extracts only asset price from VOut
	R1=VOut[1];	//Extracts only R1 from VOut
	R2=VOut[2];	//Extracts only R2 from VOut

	//EvRiskProf:: evaluates risk profile by applying beta coefficents to valuations paths
	//			   outputs to file
	double EPE = EvRiskProf(Poly, dt);

	//Record duration
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate risk profile: "<<duration<<"s"<<endl;

	cout<<"EPE = "<<EPE<<endl;

	return Fl::run();

	// //Wait for user to manually close program
	// cout<<"Press any key to exit...";
	// cin.ignore();
	// cin.get();

}

double EvRiskProf(double Poly, double dt){
	//Declare Variables necessary for valutaion
	size_t Path_Length = V[0].size();
	size_t Num_Path = V.size();
	size_t Num_Buckets = XMain[0].size();

	//Initiate Counting varaibales
	size_t Tv = 0;
	int Interval = 12;
	double tmp, tmpR1, tmpR2;

	//Out Valuation paths (Optional)
	// if(OUT == 1){
	// 	string CVANA = "tmp/VPath.dat";
	// 	ofstream CVAA;
	// 	CVAA.open(CVANA.c_str());

	// 	for(size_t i = 0 ; i < V[0].size() ; i++){
	// 		for(size_t j = 0 ; j < V.size() ; j++){
	// 			CVAA<<V[j][i]<<"\t";
	// 		}
	// 		CVAA<<endl;
	// 	}
	// }


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
					Val_Path.push_back(CoefApp(tmpBeta,tmp,tmpR1,tmpR2, Poly));
					BF<<setprecision(4)<<Tv+1<<"\t"<<tmp<<"\t"<<CoefApp(tmpBeta,tmp,tmpR1,tmpR2,Poly)<<endl;
					k = Num_Buckets;
					tmpBeta.clear();
				}
				else{}
			}
		}
		Val_T.push_back(Val_Path);
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

	vector <double> EPEmean;

	for( size_t i = 0 ; i < Val_T.size() ; i++ ){
		for( size_t j = 0 ; j < Val_T[0].size() ; j++ ){
			CVA<<i+1<<"\t"<<Val_T[i][j]<<endl;
		}
		Mean<<i+1<<"\t"<<mean(Val_T[i])<<endl;
		EPEmean.push_back(mean(Val_T[i]));
		PMean<<i+1<<"\t"<<mean(X_percent(Val_T[i],0.95))<<endl;
	}

	return mean(EPEmean);

}


double parab(vector<double> Coef, double X){
	double Y = 0;
	for(size_t i =0 ; i<Coef.size()-2 ; i++){
		Y+=pow(X,i)*Coef[i];
	}
	if(Y<0){Y=0;}
	return Y;
}

double CoefApp(vector<double> Coef, double S, double R1, double R2, double Poly){
	double Y = 0;
	for(int i =0 ; i<=Poly ; i++){
		Y+=pow(S,i)*Coef[i];
	}
	Y+=R1*Coef[Poly+1]+R2*Coef[Poly+2];
	if(Y<0){Y=0;}
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

