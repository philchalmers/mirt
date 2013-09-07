#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

RcppExport SEXP traceLinePts(SEXP Rpar, SEXP RTheta, SEXP RasMatrix, SEXP Rot) 
{
    BEGIN_RCPP

	NumericVector par(Rpar);
    NumericVector ot(Rot);
    IntegerVector asMatrix(RasMatrix);
    NumericMatrix Theta(RTheta);
    
    const int len = par.length();
    NumericVector a(Theta.ncol());
    const double u = par(len-1);
    const double g = par(len-2);
	const double d = par(len-3);
    for(int i = 0; i < Theta.ncol(); i++)
        a(i) = par(i);    
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
    const int USEOT = ot.length() > 1;
	NumericVector P(nquad);
    NumericVector Q(nquad);
	
	int i, j;
	NumericVector z(nquad);		

	//compute item trace vector
	for (j = 0; j <	nquad; j++){
		for (i = 0; i <	nfact; i++){		
			z(j) += a(i) * Theta(j,i); 
		}
		z(j) += d;
	}	
    if(USEOT){
        for (j = 0; j < nquad; j++)
            z(j) += ot(j);
    }
	for (i = 0; i < nquad; i++){ 
		P(i) = g + (u - g) * (1.0)/(1.0 + exp((-1.0)*z(i)));		        
        if(P(i) < 1e-20) P(i) = 1e-20;
        if((1.0 - P(i)) < 1e-20) P(i) = 1.0 - 1e-20;        
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
RcppExport SEXP gradedTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritemexp, SEXP Rot) 
{
    BEGIN_RCPP

    int i,j;
    NumericVector par(Rpar);	
    NumericVector ot(Rot);
	NumericMatrix Theta(RTheta);
	IntegerVector itemexp(Ritemexp);
    NumericVector a(Theta.ncol());
    for(i = 0; i < Theta.ncol(); i++)
        a(i) = par(i);
    const int ncat = par.length() - Theta.ncol();
    NumericVector d(ncat);        
    for(i = Theta.ncol(); i < par.length(); i++)
        d(i - Theta.ncol()) = par(i);
    const double nullzero = 0.0, nullone = 1.0;
    const int nquad = Theta.nrow();
	NumericMatrix Pk(nquad, ncat + 2);
	NumericMatrix P(nquad, ncat + 1);

	for(i = 0; i < nquad; i++)
        Pk(i,0) = 1.0;
    for(i = 0; i < ncat; i++)
        Pk(_,i+1) = itemTrace(a, &d(i), Theta, &nullzero, &nullone, ot); 
    if(itemexp(0)){
        for(i = (Pk.ncol()-2); i >= 0; i--)
            P(_,i) = Pk(_,i) - Pk(_,i+1);
        for(i = 0; i < P.nrow(); i++){
            for(j = 0; j < P.ncol(); j++){
                if(P(i,j) < 1e-20) P(i,j) = 1e-20;
                if((1.0 - P(i,j)) < 1e-20) P(i,j) = 1.0 - 1e-20;        
            }
        }
        return(P);
    }

    return(Pk);
	END_RCPP
}


RcppExport SEXP nominalTraceLinePts(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta, 
    SEXP RreturnNum, SEXP Rot) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector ak(Rak);
	NumericVector d(Rd);	
    NumericVector ot(Rot);
	NumericMatrix Theta(RTheta);
	IntegerVector returnNum(RreturnNum);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	const int ncat = d.length();
    const int USEOT = ot.length() > 1;
	int i,j;

	NumericMatrix Num(nquad, ncat);
	NumericMatrix P(nquad, ncat);
	NumericVector Den(nquad);
	NumericVector innerprod(nquad);

	for(i = 0; i < nquad; i++)
	    for(j = 0; j < nfact; j++)
	        innerprod(i) += Theta(i,j) * a(j);
    if(USEOT){
        for(i = 0; i < nquad; i++){
            for(j = 0; j < ncat; j++){
    	        Num(i,j) = exp(ak(j) * innerprod(i) + d(j) + ot(i));
                Den(i) += Num(i,j);
            }        
        }
    } else {
    	for(i = 0; i < nquad; i++){
    	    for(j = 0; j < ncat; j++){
    	        Num(i,j) = exp(ak(j) * innerprod(i) + d(j));
                Den(i) += Num(i,j);
            }        
        }
    }
    if(returnNum(0)) return(Num);
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++){
	        P(i,j) = Num(i,j) / Den(i);
            if(P(i,j) < 1e-20) P(i,j) = 1e-20;
            if((1.0 - P(i,j)) < 1e-20) P(i,j) = 1.0 - 1e-20;        
        }
    }

    return(P);
	END_RCPP
}
