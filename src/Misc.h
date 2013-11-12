#include<Rcpp.h>
using namespace Rcpp;

#ifndef _MISC_H
#define _MISC_H

NumericMatrix polyOuter(NumericMatrix, NumericVector,
	NumericVector, NumericVector, NumericVector, 
	NumericVector, NumericVector);

NumericVector itemTrace(NumericVector, const double *, 
        NumericMatrix, const double *, const double *, NumericVector);
        
double logit(const double *);

double antilogit(const double *);

SEXP vec2mat(std::vector<double>, const int *, const int *);

const double ABS_MAX_Z = 30;

#endif 

