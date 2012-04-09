#include"Misc.h"

SEXP logLik(SEXP Rlambdas, SEXP Rzetas, SEXP Rguess, SEXP Rtheta0,
	SEXP Rfulldata, SEXP Ritemloc, SEXP RK, SEXP RJ, SEXP RN, SEXP Rnfact,
	SEXP RestComp){

	SEXP Rreturn;
	unsigned int i, j, k, m, nfact, J, N, Ksums = 0, max = 2;
	int *itemloc,*K,*Pfulldata,*estComp;
	double *Preturn,*Plambdas,*zetas,*guess,*Ptheta0;
		
	PROTECT(Rlambdas = AS_NUMERIC(Rlambdas));
	PROTECT(Rzetas = AS_NUMERIC(Rzetas));
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(Rtheta0 = AS_NUMERIC(Rtheta0));	
	PROTECT(Rfulldata = AS_INTEGER(Rfulldata));
	PROTECT(Ritemloc = AS_INTEGER(Ritemloc));
	PROTECT(RK = AS_INTEGER(RK));	
	PROTECT(RJ = AS_INTEGER(RJ));	
	PROTECT(RN = AS_INTEGER(RN));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RestComp = AS_INTEGER(RestComp));
	Plambdas = NUMERIC_POINTER(Rlambdas);
	zetas = NUMERIC_POINTER(Rzetas);
	guess = NUMERIC_POINTER(Rguess);
	Ptheta0 = NUMERIC_POINTER(Rtheta0);	
	Pfulldata = INTEGER_POINTER(Rfulldata);
	itemloc = INTEGER_POINTER(Ritemloc);	
	K = INTEGER_POINTER(RK);
	estComp = INTEGER_POINTER(RestComp);
	J = NUMERIC_VALUE(RJ);
	N = NUMERIC_VALUE(RN);
	nfact = NUMERIC_VALUE(Rnfact);	
	for(i = 0; i < J; i++){
		Ksums += K[i]; 
		if(max < K[i]) max = K[i];
	}
	
	PROTECT(Rreturn = NEW_NUMERIC(N));
	Preturn = NUMERIC_POINTER(Rreturn);
	double a[nfact], d[max], g, lambdas[J][nfact], irt0[N],
		Plong_0[N * max];
	unsigned int loc = 0, location[J], tmpcount = 0;
	
	k = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < J; j++){
			lambdas[j][i] = Plambdas[k];
			k++;
		}
	}		
	for(i = 0; i < J; i++)		
		location[i] = itemloc[i] * N;	
	for(i = 0; i < N; i++)
		irt0[i] = 0.0;		
		
	for(unsigned int item = 0; item < J; item++){		
		if(estComp[item]){
			double dnew[nfact];
			tmpcount = 0;
			for(i = 0; i < nfact; i++){
				if(a[i] != 0.0){
					dnew[i] = zetas[i + loc];
					tmpcount += 1;
				} else dnew[i] = 200;
			}
			ProbComp(Plong_0, &k, &N, &nfact, Ptheta0, a, dnew, &g);
			m = 0;
			for(j = 0; j < k; j++){
				for(i = 0; i < N; i++){				
					if(Pfulldata[m + location[item]])
						irt0[i] += log(Plong_0[m]);													
					m++;
				}
			}	
			loc += tmpcount; 
		} else {
			k = K[item];
			for(i = 0; i < nfact; i++)
				a[i] = lambdas[item][i];		
			for(i = 0; i < (k-1); i++) 
				d[i] = zetas[i + loc];
			g = guess[item];		
			Prob(Plong_0, &k, &N, &nfact, Ptheta0, a, d, &g);		
			m = 0;
			for(j = 0; j < k; j++){
				for(i = 0; i < N; i++){				
					if(Pfulldata[m + location[item]])
						irt0[i] += log(Plong_0[m]);													
					m++;
				}
			}		
			loc += k - 1;
		}
	}	
	for(i = 0; i < N; i++)				 		
		Preturn[i] = exp(irt0[i]);
		
	UNPROTECT(12);	
	return(Rreturn);	
}


