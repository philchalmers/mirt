#include"Misc.h"

RcppExport SEXP drawThetas(SEXP Runif, SEXP Rden0, SEXP Rden1, SEXP Rlambdas, SEXP Rzetas, 
	SEXP Rguess, SEXP Rupper, SEXP Rtheta0, SEXP Rtheta1, SEXP Rfulldata, SEXP Ritemloc, SEXP RestComp)
{
    BEGIN_RCPP
	NumericVector unif(Runif);
	NumericVector den0(Rden0);
	NumericVector den1(Rden1);
	NumericMatrix lambdas(Rlambdas);
	List zetaslist(Rzetas);
	NumericVector guess(Rguess);
	NumericVector upper(Rupper);
	NumericMatrix theta0(Rtheta0);
	NumericMatrix theta1(Rtheta1);
	IntegerMatrix fulldata(Rfulldata);
	IntegerVector itemloc(Ritemloc);
	IntegerVector estComp(RestComp);

	int i, j, J, N, nfact, nzetas, istart;
	double g, u;
	J = lambdas.nrow(); 
	N = fulldata.nrow();
	nfact = lambdas.ncol();
	NumericVector accept(N), cdloglik(1);
	cdloglik.fill(0.0);	
	NumericVector zetas, a(nfact), irt0(N), irt1(N);
	NumericMatrix P_0, P_1;

	for(i = 0; i < N; i++){
		irt0(i) = 0.0;
		irt1(i) = 0.0;
	}
	//loop over items to gather log-likelihoods
	for(int item = 0; item < J; item++){
	    zetas = zetaslist[item];
		nzetas = zetas.length();
		istart = itemloc(item);
		for(i = 0; i < nfact; i++)
			a(i) = lambdas(item,i);
		g = guess(item);
		u = upper(item);

		//part comp items
		if(estComp(item)){			
			P_0 = ProbComp(theta0, a, zetas, &g, &u);			
			P_1 = ProbComp(theta1, a, zetas, &g, &u);			
			for(j = 0; j < 2; j++){
				for(i = 0; i < N; i++){				
					if(fulldata(i,j + istart)){
						irt0(i) += log(P_0(i,j));
						irt1(i) += log(P_1(i,j));
					}													
				}
			}	
		} else { //comp items
			P_0 = Prob(theta0, a, zetas, &g, &u);			
			P_1 = Prob(theta1, a, zetas, &g, &u);			
			for(j = 0; j <= nzetas; j++){
				for(i = 0; i < N; i++){				
					if(fulldata(i,j + istart)){
						irt0(i) += log(P_0(i,j));
						irt1(i) += log(P_1(i,j));
					}				
				}
			}		
		}
	}	
	for(i = 0; i < N; i++){		
		irt0(i) += log(den0(i));
		irt1(i) += log(den1(i));
		accept(i) = irt1(i) - irt0(i);		
		if(accept(i) > 0.0) accept(i) = 0.0;
		if(unif(i) < exp(accept(i))) accept(i) = 1.0;
			else accept(i) = 0.0;
	}	
	for(i = 0; i < N; i++){		
		if(accept(i)) cdloglik(0) += irt1(i);
		    else cdloglik(0) += irt0(i);
	}
	List ret;
	ret["accept"] = accept;
	ret["cdloglik"] = cdloglik;
	return(ret);
	END_RCPP
}


