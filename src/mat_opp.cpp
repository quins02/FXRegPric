#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

using namespace std;

vector< vector <vector<double> > > ExpPolyReg (vector <vector<double> > Mat, int Order, bool Constant){
	int R = Mat.size();
	int C = Mat[0].size();

	vector< vector <vector<double> > > M3;

	vector<double> V1(C,1.0);
	vector< vector<double> > V;
	for(int i = 0; i < R; i++ ){
		V.push_back(V1);
	}

	if (Constant==1){
		M3.push_back(V);
	}

	M3.push_back(Mat);

	int k=2,i=0,j=0;

	while(k<=Order){
		while(i<R){
			while(j<C){
				V[i][j]=pow((double)Mat[i][j],k);
				j++;
			}
			j=0;
			i++;
		}
		i=0;
		k++;
		M3.push_back(V);
	}

	return M3;
}

vector< vector <vector<double> > > ExpHypReg (vector <vector<double> > Mat, int Order, bool Constant){
	int R = Mat.size();
	int C = Mat[0].size();

	//double N,D;

	vector< vector <vector<double> > > M3;

	vector<double> V1(C,1.0);
	vector< vector<double> > V;
	for(int i = 0; i < R; i++ ){
		V.push_back(V1);
	}

	if (Constant==1){
		M3.push_back(V);
	}

	M3.push_back(Mat);

	int k=2,i=0,j=0;

	while(k<=Order){
		while(i<R){
			while(j<C){
				V[i][j]=1/pow((double)Mat[i][j],k);
				j++;
			}
			j=0;
			i++;
		}
		i=0;
		k++;
		M3.push_back(V);
	}

	return M3;
}




vector <vector<double> > transpose( std::vector <std::vector<double> > Mat ){
	int R = Mat.size();
	int C = Mat[0].size();

	vector<double> V;
	vector< vector<double> > M;

	for(int i = 0; i < C; i++){
		for(int j = 0; j < R; j++){
			V.push_back(Mat[j][i]);
		}
		M.push_back(V);
		V.clear();
	}
	// cout<<M.size()<<endl<<M[0].size()<<endl;

	// for(int i = 0; i < M.size(); i++){
	// 	for(int j = 0; j < M[0].size(); j++){
	// 		cout<<M[i][j]<<"\t";
	// 	}
	// 	cout<<endl;
	// }

	return M;

}

vector <vector<double> > MatMult( std::vector <std::vector<double> > A, std::vector <std::vector<double> > B ){
	int n = A.size();
	int m = A[0].size();
	int p = B[0].size();

	vector <vector<double> > AB(n,vector<double>(p,0) );

	for (int j = 0; j < p; ++j) {
		for (int k = 0; k < m; ++k) {
			for (int i = 0; i < n; ++i) {
				AB[i][j] += A[i][k]*B[k][j];
			}
		}
	}

	// for(int i = 0; i < AB.size(); i++){
	// 	for(int j = 0; j < AB[0].size(); j++){
	// 		cout<<AB[i][j]<<"\t";
	// 	}
	// 	cout<<endl;
	// }

	return AB;

}

vector <vector<double> > MatInv( std::vector <std::vector<double> > A){
	// Define the dimension n of the matrix
	// and the signum s (for LU decomposition)
	int n = A.size();
	int s,i,j;
	vector <double> tmp;
	vector <vector<double> > Inv;

	// Define all the used matrices
	gsl_matrix * m = gsl_matrix_alloc (n, n);
	gsl_matrix * inverse = gsl_matrix_alloc (n, n);
	gsl_permutation * perm = gsl_permutation_alloc (n);

	// Fill the matrix m
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			gsl_matrix_set(m,i,j,A[i][j]);
		}
	}

	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (m, perm, &s);

	// Invert the matrix m
	gsl_linalg_LU_invert (m, perm, inverse);

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			tmp.push_back(gsl_matrix_get(inverse,i,j));
		}
		Inv.push_back(tmp);
		tmp.clear();
	}

	gsl_matrix_free (m);
	gsl_matrix_free (inverse);

	return Inv;
}

vector<double> ExtCol(vector< vector<double> > v, int Col){
	int R = v.size();

	vector <double> vec;

	for(int i = 0 ; i < R ; i++){
		vec.push_back(v[i][Col]);
	}

	return vec;
}

vector< vector <double> > File2Mat(string Name){
	vector <double> Temp;
	vector< vector <double> > InMAT;
	string line;
	double V;	

	ifstream data (Name.c_str());
	
	while(getline(data,line)){
		istringstream Lin(line);
		while(Lin){
			if( Lin >> V ){ Temp.push_back(V); }
		}
		InMAT.push_back(Temp);
		Temp.clear();
	}	

	return InMAT;
}

vector< vector < double > > GenIdent(int S){
	gsl_matrix * m = gsl_matrix_alloc(S,S);
	
	gsl_matrix_set_identity(m);

	vector <double> IN(S,0);
	vector < vector < double > > v(S,IN);

	for(int i = 0 ; i<S ; i++){
		for(int j = 0 ; j<S ; j++){
			v[i][j] = gsl_matrix_get(m,i,j);
		}
	}

	return v;
}
	
vector <vector <double> > MatAdd(vector <vector <double> > A, vector <vector <double> >B){
	vector <double> IN(A[0].size(),0);
	vector <vector <double> > Out(A.size(),IN);
	for(size_t i = 0 ; i < Out.size() ; i++){
		for(size_t j = 0 ; j < Out[0].size() ; j++){
			Out[i][j]=A[i][j]+B[i][j];
		}
	}
	return Out;
}

vector <vector <double> > MatConstMult(vector <vector <double> > A, double alpha){
		vector <double> IN(A[0].size(),0);
	vector <vector <double> > Out(A.size(),IN);
	for(size_t i = 0 ; i < Out.size() ; i++){
		for(size_t j = 0 ; j < Out[0].size() ; j++){
			Out[i][j]=A[i][j]*alpha;
		}
	}
	return Out;
}
