#include<Rcpp.h>
using namespace Rcpp;

#ifndef _MISC_H
#define _MISC_H

NumericMatrix polyOuter(NumericMatrix, NumericVector,
	NumericVector, NumericVector, NumericVector, 
	NumericVector, NumericVector);

NumericVector itemTrace(NumericVector, const double *, 
        NumericMatrix, const double *, const double *, NumericVector);

#endif 

