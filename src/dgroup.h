#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

#ifndef _DGROUP_H
#define _DGROUP_H

void matrixMult(double *c, const double *a, const double *b, 
	const unsigned int *dim);

void matrixMult4(double *e, const double *a, const double *b,
	const double *c, const double *d, const unsigned int *dim);

double tr(double *a, const unsigned int *dim);

void matrixSub(double *c, const double *a,const  double *b, 
	const unsigned int *dim);

void outer(double *c, const double *a, const double *b, 
	const  unsigned int *dim);

double inner(double *a, const double *b, const double *c, 
	const unsigned int *dim);

void symMat(double *dsig, const unsigned int *nfact);

#endif

