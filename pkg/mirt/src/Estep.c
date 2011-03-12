#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX,  
	SEXP Rnfact, SEXP Rr) {
	
	SEXP list,list_names,Rr1,Rr0,Rexpected;		
	double *itemtracev,*prior,expd;
	int *X,*nfact,*r,nquad,nitems,npat,i,j,k;	
	
	//Make pointers and protect variables	
	PROTECT(Ritemtrace = AS_NUMERIC(Ritemtrace));	
	PROTECT(Rprior = AS_NUMERIC(Rprior));	
	itemtracev = NUMERIC_POINTER(Ritemtrace);
	prior = NUMERIC_POINTER(Rprior);	
	
	PROTECT(RX = AS_INTEGER(RX));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));	
	PROTECT(Rr = AS_INTEGER(Rr));
	X = INTEGER_POINTER(RX);
	nfact = INTEGER_POINTER(Rnfact);		
	r = INTEGER_POINTER(Rr);
	nquad = LENGTH(Rprior);
	nitems = LENGTH(Ritemtrace) / nquad;
	npat = LENGTH(Rr);
	
	//declare dependent arrays and initialize	
	double posterior[nquad],expected[npat],
		itemtrace[nitems][nquad],r1[nitems][nquad],
		r0[nitems][nquad];	
	int data[npat][nitems];	
	
	for	(j = 0; j < nitems; j++)
		for (i = 0; i < npat; i++)
			data[i][j] = X[i + j*npat];	
	for	(j = 0; j < nquad; j++){
		for (i = 0; i < nitems; i++){
			itemtrace[i][j] = itemtracev[i + j*nitems];
			r1[i][j] = r0[i][j] = 0;
		}					
	}
	
	// Begin main function body here				
	for (int pat = 0; pat < npat; pat++){		
		for (k = 0; k < nquad; k++)
			posterior[k] = prior[k];
			
		for (int item = 0; item < nitems; item++){
			if (data[pat][item]) {
				for (k = 0; k < nquad; k++)
					posterior[k] = posterior[k]*itemtrace[item][k];
			} else {
				for (k = 0; k < nquad; k++)
					posterior[k] = posterior[k]*(1 - itemtrace[item][k]);
			}			
		}
		
		expd = 0;
		for (i = 0; i < nquad; i++)
			expd += posterior[i];		
		expected[pat]	= expd;		
		
		for (i = 0; i < nquad; i++)
			posterior[i] = r[pat]*posterior[i]/expd;	
			
		for (int item = 0; item < nitems; item++){
			if (data[pat][item]) {
				for (k = 0; k < nquad; k++)
					r1[item][k] += posterior[k];
			} else {
				for (k = 0; k < nquad; k++)
					r0[item][k] += posterior[k];
			}			
		}
	}	//end main 		
	
	//set R objects used for return	
	PROTECT(Rr1 = allocMatrix(REALSXP,nitems,nquad));	
	PROTECT(Rr0 = allocMatrix(REALSXP,nitems,nquad));
	PROTECT(Rexpected = allocVector(REALSXP,npat));	
	for (i = 0; i < npat; i++)
		REAL(Rexpected)[i] = expected[i];
		
	for (j = 0; j < nquad; j++){
		for (i = 0; i < nitems; i++){
			REAL(Rr1)[i + j + j*(nitems-1)] = r1[i][j];
			REAL(Rr0)[i + j + j*(nitems-1)] = r0[i][j];
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
		
	UNPROTECT(10);	
	return(list);
}


//Estep for bfactor
SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, 
    SEXP RX, SEXP Rnfact, SEXP Rr, SEXP Rsitems) {
	
	SEXP list,list_names,Rr1,Rr0,Rexpected;		
	double *itemtracev,*prior,*sitems;
	int *X,*nfact,*r,nquad,npquad,nitems,npat,i,j,k,sfact;	
	
	//Make pointers and protect variables	
	PROTECT(Ritemtrace = AS_NUMERIC(Ritemtrace));	
	PROTECT(Rprior = AS_NUMERIC(Rprior));		
	PROTECT(Rsitems = AS_NUMERIC(Rsitems));
	itemtracev = NUMERIC_POINTER(Ritemtrace);
	prior = NUMERIC_POINTER(Rprior);		
	sitems = NUMERIC_POINTER(Rsitems);
	
	PROTECT(RX = AS_INTEGER(RX));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));	
	PROTECT(Rr = AS_INTEGER(Rr));
	X = INTEGER_POINTER(RX);
	nfact = INTEGER_POINTER(Rnfact);		
	r = INTEGER_POINTER(Rr);
	npquad = LENGTH(Rprior);
	nquad = (int)pow(npquad,2.0);
	nitems = LENGTH(Ritemtrace) / nquad;
	npat = LENGTH(Rr);
	sfact = *nfact - 1;	
	
	//declare dependent arrays and initialize	
	double likelihoods[sfact][nquad],L[npquad][npquad],expected[npat],
		itemtrace[nitems][nquad],r1[nitems*sfact][nquad],tempsum[npquad],
		r0[nitems*sfact][nquad],Plk[sfact][npquad],Elk[sfact][npquad],
		Pls[npquad],sitemsfull[sfact][nitems],posterior[sfact][nquad];	
	int data[npat][nitems],fact;	
	
	for	(j = 0; j < nitems; j++)
		for (i = 0; i < npat; i++)
			data[i][j] = X[i + j*npat];
	k=0;			
	for	(j = 0; j < nquad; j++){
		for (i = 0; i < nitems; i++){
			itemtrace[i][j] = itemtracev[k];
			k++;
		}
	}					
	for	(j = 0; j < nquad; j++)
		for (i = 0; i < nitems*sfact; i++)
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
				if (data[pat][item]){	
					if(sitemsfull[fact][item])
					  for (k = 0; k < nquad; k++)
					  	  likelihoods[fact][k] = likelihoods[fact][k]*itemtrace[item][k];					
				} else {
					if(sitemsfull[fact][item])
					  for (k = 0; k < nquad; k++)
						likelihoods[fact][k] = likelihoods[fact][k]*(1.0 - itemtrace[item][k]);
				}			
			}
		}			
		for(fact = 0; fact < sfact; fact++){
			k=0;
			for (j = 0; j < npquad; j++){
				tempsum[j] = 0.0;
			  for (i = 0; i < npquad; i++){
			  	L[i][j] = likelihoods[fact][k];
			  	k++;
			  }
			}
			for (j = 0; j < npquad; j++)				
			  for (i = 0; i < npquad; i++)
			  	L[i][j] = L[i][j]*prior[j];
			for (j = 0; j < npquad; j++)				
			  for (i = 0; i < npquad; i++)
			    tempsum[j] += L[j][i];
			for (i = 0; i < npquad; i++)
			  Plk[fact][i] = tempsum[i];    			
		}				
		expected[pat] = 0.0;
		for (i = 0; i < npquad; i++){
		  Pls[i] = 1.0; 		  		
			for(fact = 0; fact < sfact; fact++)
			  Pls[i] = Pls[i] * Plk[fact][i];			
			expected[pat] += Pls[i] * prior[i];  
		}				
		for(fact = 0; fact < sfact; fact++)
		  for (i = 0; i < npquad; i++)
		  	Elk[fact][i] = Pls[i] / Plk[fact][i];		  	
		for(fact = 0; fact < sfact; fact++){
		  for (i = 0; i < nquad; i++){  			  	
		    posterior[fact][i] = likelihoods[fact][i]*r[i]*Elk[fact][i % npquad] 
		      / expected[pat];
		  }
	  }		    	 	
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
		
	UNPROTECT(11);	
	return(list);
}
