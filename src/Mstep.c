
static double arraysum(const double *A1, const int *length)
{  
	double Sum = 0.0;  
	for(unsigned int j = 0; j < *length; j++)		  
		Sum += A1[j];
	return (Sum);
}

static void arrayprod2(double *Prod, const double *A1, 
  const double *A2, const int *length)
{ 	
  for(unsigned int i = 0; i < *length; i++)
    Prod[i] = A1[i] * A2[i];
}

static void arrayprod3(double *Prod, const double *A1, 
   const double *A2, const double *A3, const int *length)
{ 	
	for(unsigned int i = 0; i < *length; i++)
		Prod[i] = A1[i] * A2[i] * A3[i];
}

static void itemtrace(double *P, const double *a, 
  const double *d, const double *PTheta, const double *g, 
  const int *nfact, const int *nquad)
{	
	double z[*nquad], Theta[*nquad][*nfact];
	unsigned int i, j, k;
	for (i = 0; i < *nquad; i++)
		z[i] = 0;
	k = 0;  
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *nquad; j++){		    
			Theta[j][i] = PTheta[k];
			k++;
		 }
	}	
	//compute item trace vector
	for (j = 0; j <	*nquad; j++){		
		for (i = 0; i <	*nfact; i++)		
			z[j] += 1.702 * a[i] * Theta[j][i];  		
		z[j] += *d * 1.702;
	}
	
	for (i = 0; i < *nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
}	

//Gradient
SEXP grad(SEXP Ra, SEXP Rd, SEXP Rr1, SEXP RN, SEXP Rguess, 
	SEXP RTheta, SEXP Rprior, SEXP Rparprior) 
{
	//Protect and create vars
	SEXP Rreturn;
	PROTECT(Ra = AS_NUMERIC(Ra));
	PROTECT(Rd = AS_NUMERIC(Rd));
	PROTECT(Rr1 = AS_NUMERIC(Rr1));	
	PROTECT(RN = AS_NUMERIC(RN));		
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(RTheta = AS_NUMERIC(RTheta));	
	PROTECT(Rprior = AS_NUMERIC(Rprior));
	PROTECT(Rparprior = AS_NUMERIC(Rparprior));	
			
	double *Pguess, *Pr1, *PN, *Pa, *Pd, *PTheta, *Preturn, 
		*Pprior, *Pparprior;
	const int nfact = LENGTH(Ra);
	const int nquad = LENGTH(RTheta) / nfact;
	unsigned int i, j, k;
	Pr1 = NUMERIC_POINTER(Rr1);
	PN = NUMERIC_POINTER(RN);	
	Pguess = NUMERIC_POINTER(Rguess);
	PTheta = NUMERIC_POINTER(RTheta);
	Pa = NUMERIC_POINTER(Ra);
	Pd = NUMERIC_POINTER(Rd);
	Pprior = NUMERIC_POINTER(Rprior);
	Pparprior = NUMERIC_POINTER(Rparprior);		
	PROTECT(Rreturn = NEW_NUMERIC(nfact + 1));			
	Preturn = NUMERIC_POINTER(Rreturn);    

	double P[nquad], PQ[nquad], DIF[nquad], Theta[nquad][nfact],
		tempTheta[nquad], tempArray[nquad];
	k = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nquad; j++){		    
		   Theta[j][i] = PTheta[k];
		   k++;
		}
	}	
	itemtrace(P, Pa, Pd, PTheta, Pguess, &nfact, &nquad);	  	  		  
	for(i = 0; i < nquad; i++){      
		PQ[i] = P[i] * (1.0 - P[i]) * Pprior[i];
		DIF[i] = ((Pr1[i] / PN[i]) - P[i]) * Pprior[i];	  	
	}	
	//load gradient		  
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nquad; j++)
			tempTheta[j] = Theta[j][i]; 
		arrayprod3(tempArray, PN, DIF, tempTheta, &nquad);
		Preturn[i] = arraysum(tempArray, &nquad);		
	}	
	arrayprod2(tempArray, PN, DIF, &nquad);
	Preturn[nfact] = arraysum(tempArray, &nquad);
	
	//priors
	if(Pparprior[0] > 1.0){	  
		double c, d2 = 1.0;
		for(i = 0; i < nfact; i++)
			d2 += Pa[i] * Pa[i];	  	  
		c = 2.0 * (Pparprior[0] - 1) / d2;
		for(i = 0; i < nfact; i++)
			Preturn[i] -= c * Pa[i];
	}
	if(Pparprior[2] > 0.0){
		double normprior;
		normprior = dnorm(Pparprior[1],Pparprior[1], Pparprior[2],0) -
			dnorm(*Pd, Pparprior[1], Pparprior[2], 0);
		if(*Pd < 0.0) 
			Preturn[nfact] += 2*normprior; 
		else 
			Preturn[nfact] -=  2*normprior;
	}

	for(i = 0; i <= nfact; i++)
		Preturn[i] = (-1) * Preturn[i];
    
	UNPROTECT(9);		
	return(Rreturn);	
}

//Log-likelihood
SEXP loglik(SEXP Ra, SEXP Rd, SEXP Rr1, SEXP RN, SEXP Rguess, 
	SEXP RTheta, SEXP Rparprior) 
{
	//Proctect and create vars
	SEXP Rreturn;
	PROTECT(Ra = AS_NUMERIC(Ra));
	PROTECT(Rd = AS_NUMERIC(Rd));
	PROTECT(Rr1 = AS_NUMERIC(Rr1));	
	PROTECT(RN = AS_NUMERIC(RN));		
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(RTheta = AS_NUMERIC(RTheta));	
	PROTECT(Rparprior = AS_NUMERIC(Rparprior));	
			
	double *Pguess, *Pr1, *PN, *Pa, *Pd, *PTheta, *Preturn, 
	  *Pparprior;
	const int nfact = LENGTH(Ra);
	const int nquad = LENGTH(RTheta) / nfact;
	unsigned int i;
	Pr1 = NUMERIC_POINTER(Rr1);
	PN = NUMERIC_POINTER(RN);	
	Pguess = NUMERIC_POINTER(Rguess);
	PTheta = NUMERIC_POINTER(RTheta);
	Pa = NUMERIC_POINTER(Ra);
	Pd = NUMERIC_POINTER(Rd);
	Pparprior = NUMERIC_POINTER(Rparprior);		
	PROTECT(Rreturn = NEW_NUMERIC(1));			
	Preturn = NUMERIC_POINTER(Rreturn);    

	double P[nquad], Q[nquad], l = 0.0, sigma = 1.0;
    itemtrace(P, Pa, Pd, PTheta, Pguess, &nfact, &nquad);
	for (i = 0; i < nquad; i++)
		Q[i] = 1 - P[i];
	for(i = 0; i < nquad; i++) 
		l += Pr1[i]*log(P[i]) + (PN[i] - Pr1[i])*log(Q[i]);
    //priors
	if(Pparprior[0] > 1.0){		
		double d = 1.0, alpha[nfact];
		for (i = 0; i < nfact; i++)
			d += Pa[i]*Pa[i];
		d = pow(d,0.5);
		for (i = 0; i < nfact; i++)
			alpha[i] = Pa[i] / d;
		for (i = 0; i < nfact; i++)
			sigma -= alpha[i]*alpha[i];		
		l += log(pow(sigma,Pparprior[0] - 1.0) / beta(Pparprior[0],1.0));				
	}
	if(Pparprior[2] > 0.0)		
		l += log(dnorm(*Pd,Pparprior[1],Pparprior[2],0));     

    *Preturn = (-1.0)*l;		
	UNPROTECT(8);		
	return(Rreturn);	
}

