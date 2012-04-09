#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

#ifndef _MSTEP_H
#define _MSTEP_H

double arraysum(const double *A1, const int *length);

void arrayprod2(double *Prod, const double *A1, 
  const double *A2, const int *length);

void arrayprod3(double *Prod, const double *A1, 
   const double *A2, const double *A3, const int *length);

#endif
