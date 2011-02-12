#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP Rprior2,
    SEXP RX, SEXP Rnfact, SEXP Rr, SEXP Rsitems) {
	
	SEXP list,list_names,Rr1,Rr0,Rexpected;		
	double *itemtracev,*prior,*prior2,*sitems;
	int *X,*nfact,*r,nquad,nitems,npat,i,j,k,sfact;	
	
	//Make pointers and protect variables	
	PROTECT(Ritemtrace = AS_NUMERIC(Ritemtrace));	
	PROTECT(Rprior = AS_NUMERIC(Rprior));	
	PROTECT(Rprior2 = AS_NUMERIC(Rprior2));
	PROTECT(Rsitems = AS_NUMERIC(Rsitems));
	itemtracev = NUMERIC_POINTER(Ritemtrace);
	prior = NUMERIC_POINTER(Rprior);	
	prior2 = NUMERIC_POINTER(Rprior2);	
	sitems = NUMERIC_POINTER(Rsitems);
	
	PROTECT(RX = AS_INTEGER(RX));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));	
	PROTECT(Rr = AS_INTEGER(Rr));
	X = INTEGER_POINTER(RX);
	nfact = INTEGER_POINTER(Rnfact);		
	r = INTEGER_POINTER(Rr);
	nquad = LENGTH(Rprior);
	nitems = LENGTH(Ritemtrace) / nquad;
	npat = LENGTH(Rr);
	sfact = *nfact - 1;
	
	
	//declare dependent arrays and initialize	
	double likelihoods[sfact][nquad],expected[npat],
		itemtrace[nitems][nquad],r1[nitems*sfact][nquad],
		r0[nitems*sfact][nquad],Plk[sfact],Ek[sfact],Pls = 1.0,
		sitemsfull[sfact][nitems],posterior[sfact][nquad];	
	int data[npat][nitems],fact;	
	
	for	(j = 0; j < nitems; j++)
		for (i = 0; i < npat; i++)
			data[i][j] = X[i + j*npat];	
	for	(j = 0; j < nquad; j++)
		for (i = 0; i < nitems; i++)
			itemtrace[i][j] = itemtracev[i + j*nitems];			
	for	(j = 0; j < nquad; j++)
		for (i = 0; i < nitems*(sfact); i++)
			r1[i][j] = r0[i][j] = 0;
  k = 0; 	
	for	(j = 0; j < nitems; j++){
		for (i = 0; i < sfact; i++){
		  sitemsfull[i][j] = sitems[k];
		  k++;
		}
	}	  			
		
	// Begin main function body here				
	for (int pat = 0; pat < npat; pat++){
    for(fact = 0; fact < sfact; fact++){ 	
			for (k = 0; k < nquad; k++)
				likelihoods[fact][k] = 1;				
			for (int item = 0; item < nitems; item++){
				if (data[pat][item]) {
					for (k = 0; k < nquad; k++)
						likelihoods[fact][k] = likelihoods[fact][k]*pow(itemtrace[item][k],sitemsfull[fact][item]);
				} else {
					for (k = 0; k < nquad; k++)
						likelihoods[fact][k] = likelihoods[fact][k]*pow((1 - itemtrace[item][k]),sitemsfull[fact][item]);
				}			
			}
		}		  		  
		for(fact = 0; fact < sfact; fact++){
		  Plk[fact] = 0.0;
      for (k = 0; k < nquad; k++) 	
			  Plk[fact] += likelihoods[fact][k]*prior2[k];
		}	  
		Pls = 1.0;	  
		for(fact = 0; fact < sfact; fact++)
			Pls = Pls*Plk[fact];
    expected[pat] = Pls;			
	  for(fact = 0; fact < sfact; fact++)
			Ek[fact] = Pls/Plk[fact];		
		for(fact = 0; fact < sfact; fact++)		
			for (i = 0; i < nquad; i++)
				posterior[fact][i] = r[pat]*likelihoods[fact][i]*Ek[fact] / Pls;	
		for(fact = 0; fact < sfact; fact++){			
			for (int item = 0; item < nitems; item++){
				if (data[pat][item]) {
					for (k = 0; k < nquad; k++)
						r1[item + nitems*fact][k] += posterior[fact][k];
				} else {
					for (k = 0; k < nquad; k++)
						r0[item + nitems*fact][k] += posterior[fact][k];
				}			
			}
		}	
	}	//end main 
	
	//set R objects used for return	
	PROTECT(Rr1 = allocMatrix(REALSXP,nitems*sfact,nquad));	
	PROTECT(Rr0 = allocMatrix(REALSXP,nitems*sfact,nquad));
	PROTECT(Rexpected = allocVector(REALSXP,npat));	
	for (i = 0; i < npat; i++)
		REAL(Rexpected)[i] = expected[i];
	
	k = 0;		
	for (j = 0; j < nquad; j++){
		for (i = 0; i < nitems*sfact; i++){
			REAL(Rr1)[k] = r1[i][j];
			REAL(Rr0)[k] = r0[i][j];
			k++;
		}
	}	
		
	//set list names		
	char *names[3] = {"r1","r0","expected"};
	PROTECT(list_names = allocVector(STRSXP,3));
	for(i = 0; i < 3; i++)
		SET_STRING_ELT(list_names, i, mkChar(names[i]));
	//set list
	PROTECT(list = allocVector(VECSXP,3));
	SET_VECTOR_ELT(list, 0, Rr1);
	SET_VECTOR_ELT(list, 1, Rr0);
	SET_VECTOR_ELT(list, 2, Rexpected);		
	setAttrib(list, R_NamesSymbol, list_names); 
		
	UNPROTECT(12);	
	return(list);
}
