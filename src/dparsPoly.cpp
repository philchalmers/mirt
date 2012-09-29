#include"Misc.h"

RcppExport SEXP dparsPoly(SEXP Rprob, SEXP RThetas, SEXP RPrior, SEXP Rdat, SEXP Rnzeta) 
{		
    BEGIN_RCPP
    /* 
        Rprob = numeric matrix of probabilities
        RThetas = numeric matrix of abilities
        Rdat = integer matrix of dichotomized item responses
        nzeta = integer number of response categories
     */

	int i, j, k, nzeta, nfact, N; 
	NumericMatrix prob(Rprob);
	NumericMatrix Thetas(RThetas);
    NumericMatrix Prior(RPrior);
    NumericMatrix dat(Rdat);
    IntegerVector Pnzeta(Rnzeta);
    nzeta = Pnzeta[0];
    nfact = Thetas.ncol();
    N = Thetas.nrow();
    NumericMatrix d2L(nfact + nzeta, nfact + nzeta);
    NumericVector dL(nfact + nzeta);

	NumericVector Pk(N), Pk_1(N), Pk_p1(N), PQ_1(N), PQ(N), PQ_p1(N), 
			Pk_1Pk(N), Pk_Pkp1(N), dif1(N), dif1sq(N), dif2(N), 
			dif2sq(N), tmp1(N), tmp2(N), tmp3(N), csums(nfact);			
	NumericMatrix P(N,nzeta+2), PQfull(N,nzeta+2), mattmp(N,nfact), d2Louter;	   
	double tmp;
	IntegerVector factind(nfact);
	for(j = 0; j < (nzeta + 2); j++){
		for(i = 0; i < N; i++){
			P(i,j) = prob(i,j);
			PQfull(i,j) = prob(i,j) * (1.0 - prob(i,j));
		}
	}
	for(j = 0; j < nfact; j++)
		factind(j) = nzeta + j;
	for(j = 0; j < (nzeta + 1); j++){
		if(j < nzeta){
			for(i = 0; i < N; i++){
				Pk_1(i) = P(i,j);
				Pk(i) = P(i,j + 1);
				Pk_p1(i) = P(i,j + 2);
				PQ_1(i) = PQfull(i,j);
				PQ(i) = PQfull(i,j + 1);
				PQ_p1(i) = PQfull(i,j + 2);
				Pk_1Pk(i) = Pk_1(i) - Pk(i);
				Pk_Pkp1(i) = Pk(i) - Pk_p1(i);
				if(Pk_1Pk(i) < 1e-10) Pk_1Pk(i) = 1e-10;
				if(Pk_Pkp1(i) < 1e-10) Pk_Pkp1(i) = 1e-10;
				dif1(i) = dat(i,j) / Pk_1Pk(i);
				dif1sq(i) = dat(i,j) / (Pk_1Pk(i) * Pk_1Pk(i));
				dif2(i) = dat(i,j+1) / Pk_Pkp1(i);
				dif2sq(i) = dat(i,j+1) / (Pk_Pkp1(i) * Pk_Pkp1(i));
			}			
			tmp = 0.0;
			for(i = 0; i < N; i++)
				tmp += (-1.0) * PQ(i) * (dif1(i) - dif2(i));			
			dL(j) = tmp;			
			tmp = 0.0;
			for(i = 0; i < N; i++)
				tmp += (-1.0) * PQ(i) * PQ(i) * (dif1sq(i) + dif2sq(i)) -				
					(dif1(i) - dif2(i)) * (Pk(i) * (1.0 - Pk(i)) * (1.0 - 2.0*Pk(i)));			
			d2L(j,j) = tmp;
			if(j < (nzeta - 1)){
				tmp = 0.0;
				for(i = 0; i < N; i++)
					tmp += dif2sq(i) * PQ_p1(i) * PQ(i);
				d2L(j,j+1) = tmp;
				d2L(j+1,j) = tmp;
			}
			for(i = 0; i < N; i++){
				tmp1(i) = (-1.0) * dif2sq(i) * PQ(i) * (PQ(i) - PQ_p1(i));
				tmp2(i) = dif1sq(i) * PQ(i) * (PQ_1(i) - PQ(i));
				tmp3(i) = (dif1(i) - dif2(i)) * (Pk(i) * (1.0 - Pk(i)) * (1.0 - 2.0*Pk(i)));
			}
			for(k = 0; k < nfact; k++){
				csums(k) = 0.0;
				for(i = 0; i < N; i++){
					mattmp(i,k) = tmp1(i) * Thetas(i,k) + tmp2(i) * Thetas(i,k) - 
						tmp3(i) * Thetas(i,k);
					csums(k) += mattmp(i,k);
				}
			}
			for(i = 0; i < nfact; i++){
				d2L(j,factind(i)) = csums(i);
				d2L(factind(i),j) = csums(i);
			}			
		} else {					
			for(i = 0; i < N; i++){
				Pk_1(i) = P(i,j);
				Pk(i) = P(i,j + 1);			
				PQ_1(i) = PQfull(i,j);
				PQ(i) = PQfull(i,j + 1);			
				Pk_1Pk(i) = Pk_1(i) - Pk(i);			
				if(Pk_1Pk(i) < 1e-10) Pk_1Pk(i) = 1e-10;			
				dif1(i) = dat(i,j) / Pk_1Pk(i);
				dif1sq(i) = dat(i,j) / (Pk_1Pk(i) * Pk_1Pk(i));			
			}	
		}
		for(k = 0; k < nfact; k++){
			csums(k) = 0.0;
			for(i = 0; i < N; i++){
				mattmp(i,k) = dif1(i) * (PQ_1(i) - PQ(i)) * Thetas(i,k);
				csums(k) += mattmp(i,k);
			}
		}
		for(i = 0; i < nfact; i++)
    		dL(factind(i)) += csums(i);			
		
		d2Louter = polyOuter(Thetas, Pk, Pk_1, PQ_1, PQ, dif1sq, dif1);		
		for(k = 0; k < nfact; k++)
			for(i = 0; i < nfact; i++)
				d2L(factind(i),factind(k)) += d2Louter(i,k);				
	}

    List ret;
    ret["grad"] = dL;
    ret["hess"] = d2L;
	return(ret);
	END_RCPP
}

