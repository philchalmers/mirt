#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

#ifndef _MISC_H
#define _MISC_H

void polyOuter(double *, const double *, const double *, const double*, const double*, 
        const double *, const double *, const double *, const unsigned int *, 
        const unsigned int *);

void itemtrace(double *, const double *, const double *, const double *, const double *, 
        const unsigned int *, const unsigned int *);

void Prob(double *, const unsigned int *, const unsigned int *,	const unsigned int *, 
        const double *, const double *, const double *, const double *);

void ProbComp(double *, const unsigned int *, const unsigned int *, const unsigned int *, 
        const double *, const double *, const double *, const double *);

#endif 



