#include<Rcpp.h>
using namespace Rcpp;

#ifndef _MISC_H
#define _MISC_H

NumericMatrix polyOuter(const NumericMatrix &, const NumericVector &,
	const NumericVector &, const NumericVector &, const NumericVector &, 
	const NumericVector &, const NumericVector &);

NumericVector itemTrace(const NumericVector &, const double *, 
        const NumericMatrix &, const double *, const double *, const NumericVector &);
        
double logit(const double *);

double antilogit(const double *);

SEXP vec2mat(std::vector<double> &, const int *, const int *);

const double ABS_MAX_Z = 30;

#endif 

