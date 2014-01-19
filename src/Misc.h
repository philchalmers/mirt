#include<Rcpp.h>
using namespace Rcpp;
using std::vector;

#ifndef _MISC_H
#define _MISC_H

NumericMatrix polyOuter(const NumericMatrix &, const vector<double> &,
	const vector<double> &, const vector<double> &, const vector<double> &,
	const vector<double> &, const vector<double> &);

void itemTrace(vector<double> &, vector<double> &, const vector<double> &, const double *,
        const NumericMatrix &, const double *, const double *, const NumericVector &);

double logit(const double *);

double antilogit(const double *);

double vecsum(const vector<double> &);

SEXP vec2mat(vector<double> &, const int &, const int &);

const double ABS_MAX_Z = 30;

#endif

#ifdef _OPENMP
#include <omp.h>
#endif