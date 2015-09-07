#ifndef MAT_OPP_H
#define MAT_OPP_H

#include <vector>
#include <string>

std::vector< std::vector < std::vector<double> > > ExpPolyReg (std::vector < std::vector<double> > Mat, int Order, bool Constant);
std::vector< std::vector < std::vector<double> > > ExpHypReg (std::vector < std::vector<double> > Mat, int Order, bool Constant);
std::vector <std::vector<double> > transpose( std::vector <std::vector<double> > Mat );
std::vector <std::vector<double> > MatMult( std::vector <std::vector<double> > A, std::vector <std::vector<double> > B );
std::vector <std::vector<double> > MatInv( std::vector <std::vector<double> > A);
std::vector <double> ExtCol(std::vector<std::vector<double> > v, int Col);
std::vector <std::vector <double> > File2Mat(std::string Name);
std::vector <std::vector <double> > GenIdent(int S);
std::vector <std::vector <double> > MatAdd(std::vector <std::vector <double> > A, std::vector <std::vector <double> >B);
std::vector <std::vector <double> > MatConstMult(std::vector <std::vector <double> > A, double alpha);
#endif