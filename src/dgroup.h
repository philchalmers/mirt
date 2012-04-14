#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

#ifndef _DGROUP_H
#define _DGROUP_H

void matrixMult(double *, const double *, const double *, const unsigned int *);

void matrixMult4(double *, const double *, const double *, const double *, 
        const double *, const unsigned int *);

double tr(double *, const unsigned int *);

void matrixSub(double *, const double *, const double *, const unsigned int *);

void outer(double *, const double *, const double *, const unsigned int *);

double inner(double *, const double *, const double *, const unsigned int *);

void symMat(double *, const unsigned int *);

#endif

