#include"Misc.h"

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


