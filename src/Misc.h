#ifndef _MISC_H
#define _MISC_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using std::vector;

NumericMatrix polyOuter(const NumericMatrix &, const vector<double> &,
	const vector<double> &, const vector<double> &, const vector<double> &,
	const vector<double> &, const vector<double> &);

double logit(const double *);

double antilogit(const double *);

double vecsum(const vector<double> &);

SEXP vec2mat(vector<double> &, const int &, const int &);

const double ABS_MAX_Z = 35;

#endif
