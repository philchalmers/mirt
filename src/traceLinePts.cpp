#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

const double SQRT_DBL_MIN = sqrt(DBL_MIN);

RcppExport SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RTheta, SEXP RD) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector d(Rd);
	NumericVector g(Rg);
	NumericVector u(Ru);
	NumericVector D(RD);
	NumericMatrix Theta(RTheta);
    int nquad = Theta.nrow();
	int nfact = Theta.ncol();
	NumericVector P(nquad);
	
	int i, j;
	NumericVector z(nquad);		
	z.fill(0);

	//compute item trace vector
	for (j = 0; j <	nquad; j++){
		for (i = 0; i <	nfact; i++){		
			z(j) += D(0) * a(i) * Theta(j,i); 
		}
		z(j) += d(0) * D(0);
	}	
	for (i = 0; i < nquad; i++){ 
		P(i) = g(0) + (u(0) - g(0)) * (exp(z(i))/(1 + exp(z(i))));		
        if(P(i) < SQRT_DBL_MIN) P(i) = SQRT_DBL_MIN;
        if((1.0 - P(i)) < SQRT_DBL_MIN) P(i) = 1.0 - SQRT_DBL_MIN;        
	}
		
	return(P);

	END_RCPP
}

// graded
RcppExport SEXP gradedTraceLinePts(SEXP Ra, SEXP Rd, SEXP RTheta, SEXP RD, SEXP Ritemexp) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector d(Rd);
	NumericVector D(RD);
	NumericMatrix Theta(RTheta);
	IntegerVector itemexp(Ritemexp);
    double nullzero = 0.0, nullone = 1.0;
    int nquad = Theta.nrow();
	int nfact = Theta.ncol();
	int ncat = d.length();
	int i,j;

	NumericMatrix Pk(nquad, ncat + 2);
	NumericMatrix P(nquad, ncat + 1);

	for(i = 0; i < nquad; i++)
        Pk(i,0) = 1.0;
    for(i = 0; i < ncat; i++)
        Pk(_,i+1) = itemTrace(a, &d(i), Theta, &nullzero, &nullone, &D(0)); 
    if(itemexp(0)){
        for(i = (Pk.ncol()-2); i >= 0; i--)
            P(_,i) = Pk(_,i) - Pk(_,i+1);
        for(i = 0; i < P.nrow(); i++){
            for(j = 0; j < P.ncol(); j++){
                if(P(i,j) < SQRT_DBL_MIN) P(i,j) = SQRT_DBL_MIN;
                if((1.0 - P(i,j)) < SQRT_DBL_MIN) P(i,j) = 1.0 - SQRT_DBL_MIN;        
            }
        }
        return(P);
    }

    return(Pk);
	END_RCPP
}


RcppExport SEXP nominalTraceLinePts(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta, SEXP RD, SEXP RreturnNum) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector ak(Rak);
	NumericVector d(Rd);
	NumericVector D(RD);
	NumericMatrix Theta(RTheta);
	IntegerVector returnNum(RreturnNum);
    int nquad = Theta.nrow();
	int nfact = Theta.ncol();
	int ncat = d.length();
	int i,j;

	NumericMatrix Num(nquad, ncat);
	NumericMatrix P(nquad, ncat);
	NumericVector Den(nquad);
	NumericVector innerprod(nquad);

	for(i = 0; i < nquad; i++)
	    for(j = 0; j < nfact; j++)
	        innerprod(i) += Theta(i,j) * a(j);
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++){
	        Num(i,j) = exp(D(0) * ak(j) * innerprod(i) + D(0) * d(j));
            Den(i) += Num(i,j);
        }        
    }
    if(returnNum(0)) return(Num);
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++){
	        P(i,j) = Num(i,j) / Den(i);
            if(P(i,j) < SQRT_DBL_MIN) P(i,j) = SQRT_DBL_MIN;
            if((1.0 - P(i,j)) < SQRT_DBL_MIN) P(i,j) = 1.0 - SQRT_DBL_MIN;        
        }
    }

    return(P);
	END_RCPP
}
