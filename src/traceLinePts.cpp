#include<Rcpp.h>
using namespace Rcpp;

RcppExport SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RTheta) 
{
    BEGIN_RCPP

    /* 
        Ra = numeric vector. Item slopes 
        Rd = numeric vector. Item intercepts
        Rg = numeric scalar. Guessing parameter
        RTheta = numeric matrix. Theta values     
     */

	NumericVector a(Ra);
	NumericVector d(Rd);
	NumericVector g(Rg);
	NumericVector u(Ru);
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
			z(j) += 1.702 * a(i) * Theta(j,i); 
		}
		z(j) += d(0) * 1.702;
	}	
	for (i = 0; i < nquad; i++) 
		P(i) = g(0) + (u(0) - g(0)) * (exp(z(i))/(1 + exp(z(i))));		
		
	return(P);

	END_RCPP
}

