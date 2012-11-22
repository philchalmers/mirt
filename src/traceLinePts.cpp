#include<Rcpp.h>
using namespace Rcpp;

RcppExport SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RTheta, SEXP RD) 
{
    BEGIN_RCPP

    /* 
        Ra = numeric vector. Item slopes 
        Rd = numeric vector. Item intercepts
        Rg = numeric scalar. Guessing parameter
        RTheta = numeric matrix. Theta values     
        RD = numeric vector. Scaling parameteter
     */

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
        if(P(i) < 1e-8) P(i) = 1e-8;
        if((1.0 - P(i)) < 1e-8) P(i) = 1.0 - 1e-8;        
	}
		
	return(P);

	END_RCPP
}

