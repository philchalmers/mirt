#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

const double ABS_MAX_Z = 30;

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
    z.fill(d);

	//compute item trace vector
	for (j = 0; j <	nquad; j++){
		for (i = 0; i <	nfact; i++)		
			z(j) += a(i) * Theta(j,i); 
	}	
    if(USEOT){
        for (j = 0; j < nquad; j++)
            z(j) += ot(j);
    }
	for (i = 0; i < nquad; i++){ 
        if(z(i) > ABS_MAX_Z) z(i) = ABS_MAX_Z;
        else if(z(i) < -ABS_MAX_Z) z(i) = -ABS_MAX_Z;
		P(i) = g + (u - g) /(1.0 + exp(-z(i)));
	}
	
    if(asMatrix(0)){
        NumericMatrix ret(nquad, 2);
        ret(_, 0) = 1.0 - P;
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
                else if((1.0 - P(i,j)) < 1e-20) P(i,j) = 1.0 - 1e-20;        
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
    double z;

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
                z = ak(j) * innerprod(i) + d(j) + ot(i);
                if(z > ABS_MAX_Z) z = ABS_MAX_Z;
                else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
    	        Num(i,j) = exp(z);
                Den(i) += Num(i,j);
            }        
        }
    } else {
    	for(i = 0; i < nquad; i++){
    	    for(j = 0; j < ncat; j++){
                z = ak(j) * innerprod(i) + d(j);
                if(z > ABS_MAX_Z) z = ABS_MAX_Z;
                else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
    	        Num(i,j) = exp(z);
                Den(i) += Num(i,j);
            }        
        }
    }
    if(returnNum(0)) return(Num);
	for(i = 0; i < nquad; i++){
	    for(j = 0; j < ncat; j++)
	        P(i,j) = Num(i,j) / Den(i);
    }

    return(P);
	END_RCPP
}
