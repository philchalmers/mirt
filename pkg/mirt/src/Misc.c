#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

SEXP traceLinePts(SEXP Ra, SEXP Rd, SEXP Rg, 
	SEXP RTheta, SEXP Rnquad, SEXP Rnfact) {
	
	SEXP Rreturn;			
	int *nquad,*nfact;
	double *a,*d,*g,*Theta,*P;	
	
	//set pointers and protect
	PROTECT(Rnquad = AS_INTEGER(Rnquad));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	nquad = INTEGER_POINTER(Rnquad);
	nfact = INTEGER_POINTER(Rnfact);
	
	double z[*nquad];		
	for (int i = 0; i < *nquad; i++)
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
	for (int j = 0; j <	*nquad; j++){
		for (int i = 0; i <	*nfact; i++){		
			z[j] += 1.702 * a[i] * Theta[j + i*(*nquad)]; 
		}
		z[j] += *d * 1.702;
	}
	
	for (int i = 0; i < *nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
		
	UNPROTECT(7);	
	return(Rreturn);
}

//
SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP Rnfact, SEXP RN){
	
	SEXP Rreturn;			
	int i, j, n, nfact, N;
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
	int i, j, n, nfact, N;
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