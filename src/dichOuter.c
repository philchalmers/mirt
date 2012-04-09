#include"Misc.h"

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


