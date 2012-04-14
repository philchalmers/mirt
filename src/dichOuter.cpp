#include <Rcpp.h>
using namespace Rcpp;

// NOTE REMOVE Rnfact in R
RcppExport SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP RN){
	
	int i, j, n, nfact;
    NumericMatrix Thetas(RThetas);    
    NumericVector PQ(RPQ);
    NumericVector N(RN);
    nfact = Thetas.ncol();
	NumericMatrix ret(nfact,nfact);			

	for(n = 0; n < N[0]; n++)
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				ret(i,j) += Thetas(n,i) * Thetas(n,j) * PQ[n];
		
	return(ret);
}


