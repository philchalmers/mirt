#include"Misc.h"

RcppExport SEXP logLik(SEXP Rlambdas, SEXP Rzetas, SEXP Rguess, SEXP Rtheta0,
	SEXP Rfulldata, SEXP Ritemloc, SEXP RK,	SEXP RestComp)
{ 
    BEGIN_RCPP
	int i, j, nfact, J, N;
    NumericMatrix lambdas(Rlambdas);
    List zetas(Rzetas);
    NumericVector guess(Rguess);
    NumericMatrix theta0(Rtheta0);
    IntegerMatrix fulldata(Rfulldata);
    IntegerVector itemloc(Ritemloc); //itemloc - 1 from R
    IntegerVector K(RK);
    IntegerVector estComp(RestComp);
    NumericVector d;
    nfact = theta0.ncol();
    J = K.length();
    N = theta0.nrow();
	
	NumericVector a(nfact), irt0(N);
	double g;
	for(i = 0; i < N; i++)
		irt0(i) = 0.0;		
		
	for(int item = 0; item < J; item++){		
	    d = zetas[item];
        NumericMatrix P(N, K(item));
        for(i = 0; i < nfact; i++)
			a(i) = lambdas(item,i);
		g = guess(item);		
		if(estComp(item))
			P = ProbComp(theta0, a, d, &g); 
        else 
            P = Prob(theta0, a, d, &g);
		for(j = 0; j < K(item); j++)
			for(i = 0; i < N; i++)				
			    if(fulldata(i,j + itemloc(item)))
					irt0(i) += log(P(i,j));													
	}	

    NumericVector ret(N);
	for(i = 0; i < N; i++)				 		
		ret(i) = exp(irt0[i]);		
	return(ret);	
	END_RCPP
}


