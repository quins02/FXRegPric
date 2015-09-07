#ifndef OPT_EVAL_H
#define OPT_EVAL_H

#include <vector>

double opt_call(double S, double K, double r, double t);
double opt_put(double S, double K, double r, double t);
double opt_dig_call(double S, double K, double r, double t);
double opt_dig_put(double S, double K, double r, double t);
double barrier_call(std::vector <double> Path, double barrier, bool IO, bool UD, double K, double r, double t);
double barrier_put(std::vector <double> Path, double barrier, bool IO, bool UD, double K, double r, double t);
double PRDC(std::vector <double> FX, std::vector <double> r1, std::vector <double> r2, double N, double t);
double NORMFUNC(double x, void * params);
double NORMDIST(double x);
double call_CF(double S, double K, double R, double t, double vol);
double dig_call_CF(double S, double K, double R, double t, double vol);

#endif