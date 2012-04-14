static void Prob(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g)
{
	double Ps[*N][*k + 1], Pdif[*N][*k], p1[*N], tmp;
	unsigned int i, m = 0;
	int j;

	for(i = 0; i < *N; i++){
		Ps[i][0] = 1;
		Ps[i][*k] = 0;
	}
	for(j = 0; j < (*k - 1); j++){
		tmp = d[j];
		itemtrace(p1, a, &tmp, theta, g, nfact, N);
		for(i = 0; i < *N; i++)
			Ps[i][j + 1] = p1[i];
	}
	for(j = (*k - 1); j >= 0; j--)
		for(i = 0; i < *N; i++)
			Pdif[i][j] = Ps[i][j] - Ps[i][j + 1];				
	for(j = 0; j < *k; j++){
		for(i = 0; i < *N; i++){
			if(Pdif[i][j] < .00000001) Pdif[i][j] = .00000001;
			if(*k == 2) Pdif[i][j] = 1 - Pdif[i][j];
			P[m] = Pdif[i][j];
			m++;
		}
	}
}

static void ProbComp(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g)
{
	double Theta[*N], tmp[*N], zerog = 0.0, tmpa, tmpd;
	unsigned int i, j, onenfact = 1;	
	for(j = 0; j < *N; j++)
		tmp[j] = 1.0;
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *N; j++)
			Theta[j] = theta[j + i*(*N)];
		tmpa = a[i];
		tmpd = d[i];
		itemtrace(P, &tmpa, &tmpd, Theta, &zerog, &onenfact, N);
		for(j = 0; j < *N; j++)
			tmp[j] *= P[j];
	}
	for(j = 0; j < *N; j++){
		P[j] = *g + (1 - *g)*tmp[j];
		P[j + *N] = 1.0 - P[j];
	}

SEXP drawThetas(SEXP Runif, SEXP Rden0, SEXP Rden1, SEXP Rlambdas, SEXP Rzetas, 
	SEXP Rguess, SEXP Rtheta0, SEXP Rtheta1, SEXP Rfulldata, SEXP Ritemloc,
	SEXP RK, SEXP RJ, SEXP RN, SEXP Rnfact, SEXP RestComp){

	SEXP Rreturn;
	unsigned int i, j, k, m, nfact, J, N, Ksums = 0, max = 2;
	int *itemloc,*K,*Pfulldata,*estComp;
	double *Preturn,*Plambdas,*zetas,*guess,*Ptheta0,*Ptheta1,*unif, 
		*den0,*den1;
		
	PROTECT(Runif = AS_NUMERIC(Runif));
	PROTECT(Rden0 = AS_NUMERIC(Rden0));
	PROTECT(Rden1 = AS_NUMERIC(Rden1));
	PROTECT(Rlambdas = AS_NUMERIC(Rlambdas));
	PROTECT(Rzetas = AS_NUMERIC(Rzetas));
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(Rtheta0 = AS_NUMERIC(Rtheta0));
	PROTECT(Rtheta1 = AS_NUMERIC(Rtheta1));
	PROTECT(Rfulldata = AS_INTEGER(Rfulldata));
	PROTECT(Ritemloc = AS_INTEGER(Ritemloc));
	PROTECT(RK = AS_INTEGER(RK));	
	PROTECT(RJ = AS_INTEGER(RJ));	
	PROTECT(RN = AS_INTEGER(RN));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RestComp = AS_INTEGER(RestComp));
	unif = NUMERIC_POINTER(Runif);
	den0 = NUMERIC_POINTER(Rden0);
	den1 = NUMERIC_POINTER(Rden1);
	Plambdas = NUMERIC_POINTER(Rlambdas);
	zetas = NUMERIC_POINTER(Rzetas);
	guess = NUMERIC_POINTER(Rguess);
	Ptheta0 = NUMERIC_POINTER(Rtheta0);
	Ptheta1 = NUMERIC_POINTER(Rtheta1);
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
	
	PROTECT(Rreturn = NEW_NUMERIC(N + 1));
	Preturn = NUMERIC_POINTER(Rreturn);
	double a[nfact], d[max], g, lambdas[J][nfact], irt0[N], irt1[N], 
		accept[N], Plong_1[N * max], Plong_0[N * max], cdloglik = 0;
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
	for(i = 0; i < N; i++){
		irt0[i] = 0.0;
		irt1[i] = 0.0;
	}	
	for(unsigned int item = 0; item < J; item++){
		k = K[item];
		for(i = 0; i < nfact; i++)
			a[i] = lambdas[item][i];
		g = guess[item];
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
			ProbComp(Plong_1, &k, &N, &nfact, Ptheta1, a, dnew, &g);			
			m = 0;
			for(j = 0; j < k; j++){
				for(i = 0; i < N; i++){				
					if(Pfulldata[m + location[item]]){
						irt0[i] += log(Plong_0[m]);
						irt1[i] += log(Plong_1[m]);
					}													
					m++;
				}
			}	
			loc += tmpcount; 
		} else {					
			for(i = 0; i < (k-1); i++) 
				d[i] = zetas[i + loc];			
			Prob(Plong_0, &k, &N, &nfact, Ptheta0, a, d, &g);			
			Prob(Plong_1, &k, &N, &nfact, Ptheta1, a, d, &g);			
			m = 0;
			for(j = 0; j < k; j++){
				for(i = 0; i < N; i++){				
					if(Pfulldata[m + location[item]]){
						irt0[i] += log(Plong_0[m]);
						irt1[i] += log(Plong_1[m]);
					}				
					m++;
				}
			}		
			loc += k - 1;
		}
	}	
	for(i = 0; i < N; i++){		
		irt0[i] += log(den0[i]);
		irt1[i] += log(den1[i]);
		accept[i] = irt1[i] - irt0[i];		
		if(accept[i] > 0.0) accept[i] = 0.0;
		if(unif[i] < exp(accept[i])) accept[i] = 1.0;
			else accept[i] = 0.0;
		Preturn[i] = accept[i];
	}	
	for(i = 0; i < N; i++){		
		if(accept[i]) cdloglik += irt1[i];
		else cdloglik += irt0[i];
	}
	Preturn[N] = cdloglik;
	
	UNPROTECT(16);	
	return(Rreturn);	
}


