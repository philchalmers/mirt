#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, 
	SEXP RTheta, SEXP Rnquad, SEXP Rnfact) {
	
	SEXP Rreturn;			
	unsigned int i, j;
	int *nquad,*nfact;
	double *a,*d,*g,*Theta,*P;

	
	//set pointers and protect
	PROTECT(Rnquad = AS_INTEGER(Rnquad));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	nquad = INTEGER_POINTER(Rnquad);
	nfact = INTEGER_POINTER(Rnfact);
	
	double z[*nquad];		
	for (i = 0; i < *nquad; i++)
		z[i] = 0;		
	PROTECT(Ra = AS_NUMERIC(Ra));
	PROTECT(Rd = AS_NUMERIC(Rd));
	PROTECT(Rg = AS_NUMERIC(Rg));
	PROTECT(RTheta = AS_NUMERIC(RTheta));
	a = NUMERIC_POINTER(Ra);
	d = NUMERIC_POINTER(Rd);
	g = NUMERIC_POINTER(Rg);
	Theta = NUMERIC_POINTER(RTheta);		
	
	PROTECT(Rreturn = NEW_NUMERIC(*nquad));		
	P = NUMERIC_POINTER(Rreturn);
	
	//compute item trace vector
	for (j = 0; j <	*nquad; j++){
		for (i = 0; i <	*nfact; i++){		
			z[j] += 1.702 * a[i] * Theta[j + i*(*nquad)]; 
		}
		z[j] += *d * 1.702;
	}
	
	for (i = 0; i < *nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
		
	UNPROTECT(7);	
	return(Rreturn);
}

//
SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP Rnfact, SEXP RN){
	
	SEXP Rreturn;			
	unsigned int i, j, n, nfact, N;
	double *Preturn,*PThetas,*PQ;	
	
	//set pointers and protect
	PROTECT(RThetas = AS_NUMERIC(RThetas));	
	PROTECT(RPQ = AS_NUMERIC(RPQ));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RN = AS_INTEGER(RN));
	PThetas = NUMERIC_POINTER(RThetas);
	PQ = NUMERIC_POINTER(RPQ);
	nfact = NUMERIC_VALUE(Rnfact);
	N = NUMERIC_VALUE(RN);
	PROTECT(Rreturn = allocMatrix(REALSXP,nfact,nfact));		
	Preturn = NUMERIC_POINTER(Rreturn);

	double out[nfact][nfact], Thetas[N][nfact];
	for(i = 0; i < nfact; i++)
		for(j = 0; j < nfact; j++)
			out[i][j] = 0.0;
	for(i = 0; i < nfact; i++)
		for(n = 0; n < N; n++)
			Thetas[n][i] = PThetas[n + N*i];	
	for(n = 0; n < N; n++)
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				out[i][j] += Thetas[n][i] * Thetas[n][j] * PQ[n];
	n = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nfact; j++){
			Preturn[n] = out[j][i];	
			n++;
		}
	}
		
	UNPROTECT(5);	
	return(Rreturn);
}

SEXP polyOuter(SEXP RThetas, SEXP RPk, SEXP RPk_1, SEXP RPQ_1, 
	SEXP RPQ, SEXP Rdif1sq, SEXP Rdif1, SEXP Rnfact, SEXP RN){
	
	SEXP Rreturn;			
	unsigned int i, j, n, nfact, N;
	double *Preturn,*PThetas,*Pk,*Pk_1,*PQ,*PQ_1,*dif1sq,*dif1;		
	
	PROTECT(RThetas = AS_NUMERIC(RThetas));
	PROTECT(RPk = AS_NUMERIC(RPk));
	PROTECT(RPk_1 = AS_NUMERIC(RPk_1));
	PROTECT(RPQ = AS_NUMERIC(RPQ));
	PROTECT(RPQ_1 = AS_NUMERIC(RPQ_1));
	PROTECT(Rdif1sq = AS_NUMERIC(Rdif1sq));
	PROTECT(Rdif1 = AS_NUMERIC(Rdif1));	
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RN = AS_INTEGER(RN));
	PThetas = NUMERIC_POINTER(RThetas);
	Pk = NUMERIC_POINTER(RPk);
	Pk_1 = NUMERIC_POINTER(RPk_1);
	PQ = NUMERIC_POINTER(RPQ);
	PQ_1 = NUMERIC_POINTER(RPQ_1);
	dif1sq = NUMERIC_POINTER(Rdif1sq);
	dif1 = NUMERIC_POINTER(Rdif1);	
	nfact = NUMERIC_VALUE(Rnfact);
	N = NUMERIC_VALUE(RN);
	PROTECT(Rreturn = allocMatrix(REALSXP,nfact,nfact));		
	Preturn = NUMERIC_POINTER(Rreturn);

	double out[nfact][nfact], Thetas[N][nfact], outer[nfact][nfact],
		temp[nfact];
	for(i = 0; i < nfact; i++)
		for(j = 0; j < nfact; j++)
			out[i][j] = 0.0;
	for(i = 0; i < nfact; i++)
		for(n = 0; n < N; n++)
			Thetas[n][i] = PThetas[n + N*i];	
	for(n = 0; n < N; n++){
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				outer[i][j] = Thetas[n][i] * Thetas[n][j];
		for(i = 0; i < nfact; i++)
			temp[i] =  (PQ_1[n] * Thetas[n][i] - PQ[n] * Thetas[n][i]);
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)				
				out[i][j] += (-1) * dif1sq[n] * temp[i] * temp[j] +  
				(dif1[n] * (Pk_1[n] * (1.0 - Pk_1[n]) * (1.0 - 2.0 * Pk_1[n]) * 
				outer[i][j] - Pk[n] * (1.0 - Pk[n]) * (1.0 - 2.0 * Pk[n]) * outer[i][j]));
	}				
	n = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nfact; j++){
			Preturn[n] = out[j][i];	
			n++;
		}
	}
		
	UNPROTECT(10);	
	return(Rreturn);
}


static void itemtrace(double *P, const double *a, 
  const double *d, const double *PTheta, const double *g, 
  const int *nfact, const int *nquad)
{	
	double z[*nquad];
	unsigned int i, j, loc[*nfact];

	for (int i = 0; i < *nquad; i++)
		z[i] = 0;		
	for(i = 0; i < (*nfact); i++)
		loc[i] = i * (*nquad);
	for (j = 0; j <	*nquad; j++){		
		for (i = 0; i <	*nfact; i++){		
			z[j] += 1.702 * a[i] * PTheta[loc[i]];  		
			loc[i]++;
		}
		z[j] += *d * 1.702;
	}	
	for (i = 0; i < *nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
}

static void Prob(double *P, const int *k, const int *N, const int *nfact,
	const double *theta, const double *a, const double *d, const double *g)
{
	double Ps[*N][*k + 1], Pdif[*N][*k], p1[*N], tmp;
	int i, j;

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
	int m = 0;
	for(j = 0; j < *k; j++){
		for(i = 0; i < *N; i++){
			if(Pdif[i][j] < .00000001) Pdif[i][j] = .00000001;
			if(*k == 2) Pdif[i][j] = 1 - Pdif[i][j];
			P[m] = Pdif[i][j];
			m++;
		}
	}
}


SEXP drawThetas(SEXP Runif, SEXP Rden0, SEXP Rden1, SEXP Rlambdas, SEXP Rzetas, 
	SEXP Rguess, SEXP Rtheta0, SEXP Rtheta1, SEXP Rfulldata, SEXP Ritemloc,
	SEXP RK, SEXP RJ, SEXP RN, SEXP Rnfact){

	SEXP Rreturn;	
	int i, j, k, m, J, nfact, N,*itemloc,*K,*Pfulldata, Ksums = 0, max = 2;
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
		accept[N], Plong_0[N * k], Plong_1[N * k], cdloglik = 0;
	unsigned int loc = 0, location[J];
	
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
		for(i = 0; i < (k-1); i++) 
			d[i] = zetas[i + loc];
		g = guess[item];		
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
	
	UNPROTECT(15);	
	return(Rreturn);	
}

