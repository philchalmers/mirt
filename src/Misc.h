#include<Rcpp.h>
using namespace Rcpp;

#ifndef _MISC_H
#define _MISC_H

void polyOuter(double *, NumericMatrix, const double *,
	const double *, const double *,	const double *, 
	const double *, const double *);


void itemtrace(double *, const double *, const double *, 
        NumericMatrix, const double *);

NumericMatrix Prob(NumericMatrix, const double *,
        NumericVector, const double *);

NumericMatrix ProbComp(NumericMatrix, const double *, 
        NumericVector, const double *);

#endif 



