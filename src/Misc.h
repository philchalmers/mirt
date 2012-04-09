#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

#ifndef _MISC_H
#define _MISC_H

void polyOuter(double *d2Louter, const double *PThetas, const double *Pk,
	const double *Pk_1, const double *PQ_1,	const double *PQ, 
	const double *dif1sq, const double *dif1, const unsigned int *nfact, 
	const unsigned int *N);

void itemtrace(double *P, const double *a, 
  const double *d, const double *PTheta, const double *g, 
  const unsigned int *nfact, const unsigned int *nquad);

void Prob(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g);

void ProbComp(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g);

#endif 



