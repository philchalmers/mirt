#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

RcppExport SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RTheta, SEXP RD) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector d(Rd);
	NumericVector g(Rg);
	NumericVector u(Ru);
	NumericVector D(RD);
	NumericMatrix Theta(RTheta);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	NumericVector P(nquad);
	
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

RcppExport SEXP mcmTraceLinePts(SEXP Ra, SEXP Rak, SEXP Rd, SEXP Rt, SEXP RTheta, SEXP RD) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector ak(Rak);
	NumericVector d(Rd);
	NumericVector t(Rt);
	NumericVector D(RD);
	NumericMatrix Theta(RTheta);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	const int ncat = t.length();
	int i,j;

	NumericMatrix Num(nquad, ncat);
	NumericMatrix Num0(nquad, ncat + 1);
	NumericMatrix P(nquad, ncat);
	NumericVector Den(nquad);
	NumericVector Den0(nquad);
	NumericVector innerprod(nquad);
	NumericVector C0(nquad);

	for(i = 0; i < nquad; i++){
	    for(j = 0; j < nfact; j++)
	        innerprod(i) += Theta(i,j) * a(j);
	    Num0(i,0) = exp(D(0) * ak(0) * innerprod(i) + D(0) * d(0));
	    Den0(i) = Num0(i,0);
    }
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++){
	        Num(i,j) = exp(D(0) * ak(j+1) * innerprod(i) + D(0) * d(j+1));
	        Num0(i,j+1) = Num(i,j);
            Den(i) += Num(i,j);
            Den0(i) += Num0(i,j+1);
        }        
    }
    C0 = Num0(_,0) / Den0;
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++){
	        P(i,j) = C0(i) * t(j) + (1.0 - C0(i)) * Num(i,j) / Den(i);
            if(P(i,j) < 1e-10) P(i,j) = 1e-10;
            if((1.0 - P(i,j)) < 1e-10) P(i,j) = 1.0 - 1e-10;        
        }
    }

    return(P);
	END_RCPP
}

