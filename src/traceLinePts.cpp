#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

RcppExport SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RTheta, SEXP RD, SEXP RasMatrix) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector d(Rd);
	NumericVector g(Rg);
	NumericVector u(Ru);
	NumericVector D(RD);
    IntegerVector asMatrix(RasMatrix);
	NumericMatrix Theta(RTheta);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	NumericVector P(nquad);
    NumericVector Q(nquad);
	
	int i, j;
	NumericVector z(nquad);		

	//compute item trace vector
	for (j = 0; j <	nquad; j++){
		for (i = 0; i <	nfact; i++){		
			z(j) += D(0) * a(i) * Theta(j,i); 
		}
		z(j) += d(0) * D(0);
	}	
	for (i = 0; i < nquad; i++){ 
		P(i) = g(0) + (u(0) - g(0)) * (1.0)/(1.0 + exp((-1.0)*z(i)));		        
        if(P(i) < 1e-10) P(i) = 1e-10;
        if((1.0 - P(i)) < 1e-10) P(i) = 1.0 - 1e-10;        
        Q(i) = 1.0 - P(i);
	}
	
    if(asMatrix(0)){
        NumericMatrix ret(nquad, 2);
        ret(_, 0) = Q;
        ret(_, 1) = P;
        return(ret);
    } else return(P);

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
    const double nullzero = 0.0, nullone = 1.0;
    const int nquad = Theta.nrow();
	const int ncat = d.length();
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
                if(P(i,j) < 1e-10) P(i,j) = 1e-10;
                if((1.0 - P(i,j)) < 1e-10) P(i,j) = 1.0 - 1e-10;        
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
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	const int ncat = d.length();
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
            if(P(i,j) < 1e-10) P(i,j) = 1e-10;
            if((1.0 - P(i,j)) < 1e-10) P(i,j) = 1.0 - 1e-10;        
        }
    }

    return(P);
	END_RCPP
}
