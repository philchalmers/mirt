#include"Misc.h"

RcppExport SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP RN)
{	
    BEGIN_RCPP
	int i, j, n;
    NumericMatrix Thetas(RThetas);    
    NumericVector PQ(RPQ);
    NumericVector N(RN);
    const int nfact = Thetas.ncol();
	NumericMatrix ret(nfact,nfact);			

	for(n = 0; n < N(0); n++)
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				ret(i,j) += Thetas(n,i) * Thetas(n,j) * PQ(n);
		
	return(ret);
	END_RCPP
}

NumericMatrix polyOuter(NumericMatrix Thetas, NumericVector Pk,
	NumericVector Pk_1, NumericVector PQ_1, NumericVector PQ, 
	NumericVector dif1sq, NumericVector dif1)
{
	int i, j, n;
	const int nfact = Thetas.ncol();
	NumericMatrix d2Louter(nfact,nfact), outer(nfact,nfact);
	NumericVector temp(nfact);
	d2Louter.fill(0.0);
	
	for(n = 0; n < Thetas.nrow(); n++){
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				outer(i,j) = Thetas(n,i) * Thetas(n,j);
		for(i = 0; i < nfact; i++)
			temp(i) =  (PQ_1(n) * Thetas(n,i) - PQ(n) * Thetas(n,i));
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)				
				d2Louter(i,j) += (-1) * dif1sq(n) * temp(i) * temp(j) +  
				    (dif1(n) * (Pk_1(n) * (1.0 - Pk_1(n)) * (1.0 - 2.0 * Pk_1(n)) * 
				    outer(i,j) - Pk(n) * (1.0 - Pk(n)) * (1.0 - 2.0 * Pk(n)) * outer(i,j)));
	}
	return d2Louter;		
}

NumericVector itemTrace(NumericVector a, const double *d, 
        NumericMatrix Theta, const double *g, const double *u, const double *D)
{	
	int i, j;
    const int nquad = Theta.nrow();
	NumericVector P(nquad), z(nquad);

	for (i = 0; i <	nquad; i++){
	    z(i) = 0.0;
		for (j = 0; j <	Theta.ncol(); j++){		
			z(i) += *D * a(j) * Theta(i,j);  		
		}
		z(i) += *d * *D;
	}	
	for (i = 0; i < nquad; i++) 
		P(i) = *g + (*u - *g) * (exp(z(i))/(1 + exp(z(i))));
	
	return P;		
}

RcppExport SEXP reloadPars(SEXP Rlongpars, SEXP Rpars, SEXP Rngroups, SEXP RJ)
{    
    BEGIN_RCPP
	NumericVector longpars(Rlongpars);
    List pars(Rpars);
    NumericVector ngroups(Rngroups);
    NumericVector J(RJ);
    int i, j, g;
    int ind = 0, len;

    for(g = 0; g < ngroups[0]; g++){
        List glist = pars[g];
        for(i = 0; i < J[0]; i++){
            S4 item = glist[i];
            NumericVector p = item.slot("par");
            len = p.length();
            for(j = 0; j < len; j++)
                p(j) = longpars(ind+j);
            ind += len;
            item.slot("par") = p;
            glist[i] = item;
        }        
        pars[g] = glist;
    }
	
    return(pars);
	END_RCPP
}