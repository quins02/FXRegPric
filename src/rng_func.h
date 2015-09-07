#ifndef RNG_FUNC_H
#define RNG_FUNC_H

#include <vector>
#include <string>

std::vector <std::vector <double> > rand_gsl(int seed, int path, int T, double dt);
std::vector <std::vector <double> > CorMat(std::vector <std::vector <double> > Data);
std::vector< std::vector< std::vector <double> > > PathGen(double seed, int PATH, double T, double dt, double tmp1);

#endif