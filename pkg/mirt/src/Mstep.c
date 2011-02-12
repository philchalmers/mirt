#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>
#include <R_ext/Lapack.h>      

static double arraysum(const double *A1, const int *length)
{  
  double Sum = 0.0;  
  for(int j = 0; j < *length; j++)		  
		Sum += A1[j];
	return (Sum);
}

static void arrayprod2(double *Prod, const double *A1, 
  const double *A2, const int *length)
{ 	
  for(int i = 0; i < *length; i++)
    Prod[i] = A1[i] * A2[i];
}

static void arrayprod3(double *Prod, const double *A1, 
   const double *A2, const double *A3, const int *length)
{ 	
  for(int i = 0; i < *length; i++)
    Prod[i] = A1[i] * A2[i] * A3[i];
}

static void itemtrace(double *P, const double *a, 
  const double *d, const double *PTheta, const double *g, 
  const int *nfact, const int *nquad)
{	
	double z[*nquad], Theta[*nquad][*nfact];
	for (int i = 0; i < *nquad; i++)
		z[i] = 0;
	int k = 0;  
	for(int i = 0; i < *nfact; i++){
    for(int j = 0; j < *nquad; j++){		    
	    Theta[j][i] = PTheta[k];
	    k++;
	  }
	}	
	//compute item trace vector
	for (int j = 0; j <	*nquad; j++){		
		for (int i = 0; i <	*nfact; i++)		
			z[j] += 1.702 * a[i] * Theta[j][i];  		
		z[j] += *d * 1.702;
	}
	
	for (int i = 0; i < *nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
}	

SEXP Mstep(SEXP Rr1, SEXP RN, SEXP Rprior, SEXP Rpars,  
	SEXP Rguess, SEXP RTheta, SEXP Rparprior) 
{
	//Proctect and create vars
	SEXP Rreturn;	
	PROTECT(Rr1 = AS_NUMERIC(Rr1));	
	PROTECT(RN = AS_NUMERIC(RN));
  PROTECT(Rprior = AS_NUMERIC(Rprior));
	PROTECT(Rpars = AS_NUMERIC(Rpars));	
	PROTECT(Rguess = AS_NUMERIC(Rguess));
  PROTECT(RTheta = AS_NUMERIC(RTheta));	
	PROTECT(Rparprior = AS_NUMERIC(Rparprior));	
		
	double *guess, *Pr1, *PN, *Ppars, *PTheta, *Preturn, 
	  *Pprior, *Pparprior;
	int i, j, k;
	const int nitems = LENGTH(Rguess);
	const int nquad = LENGTH(RN) / nitems;
  const int nfact = (LENGTH(Rpars) / nitems) - 1;
  const int npars = nfact + 1; 	
  Pr1 = NUMERIC_POINTER(Rr1);
	PN = NUMERIC_POINTER(RN);
	Pprior = NUMERIC_POINTER(Rprior);
	Ppars = NUMERIC_POINTER(Rpars);
	guess = NUMERIC_POINTER(Rguess);
	PTheta = NUMERIC_POINTER(RTheta);
	Pparprior = NUMERIC_POINTER(Rparprior);
	PROTECT(Rreturn = NEW_NUMERIC(nitems * npars));			
	Preturn = NUMERIC_POINTER(Rreturn);			
	
//	//define and load arrays
	double fullr1[nitems][nquad], fullN[nitems][nquad],
	  fullpars[nitems][npars], fullparprior[nitems][3],
	  Theta[nquad][nfact];
	k = 0;  
	for(j = 0; j < nquad; j++){
    for(i = 0; i < nitems; i++){	
	    fullr1[i][j] = Pr1[k];
	    fullN[i][j] = PN[k];
	    Theta[j][i] = PTheta[k];
	    k++;
	  }
	}
	k = 0;  
	for(i = 0; i < nfact; i++){
    for(j = 0; j < nquad; j++){		    
	    Theta[j][i] = PTheta[k];
	    k++;
	  }
	}
	k = 0;
	for(j = 0; j < npars; j++){
    for(i = 0; i < nitems; i++){
      fullpars[i][j] = Ppars[k];
		  k++;
    } 	  
	}
	k = 0;
  for(j = 0; j < 3 ; j++){
    for(i = 0; i < nitems; i++){
      fullparprior[i][j] = Pparprior[k];	  
	    k++;
    }	
	}  
		
	//****************************************************
	// Main big loop
	double a[nfact], d, g, P[nquad], PQ[nquad], 
	  DIF[nquad], L[npars], temp, temparray[nquad],
	  tempTheta[nquad], r1[nquad], N[nquad], LL[npars][npars],
	  LLinv[npars][npars], correction[npars], StepLimit[npars], 
	  LastL[npars], work[npars], corSign[npars], d2, c, 
	  betaprior[nfact][nfact], normprior; 	
	int iter = 1, info = 0, ipiv[npars];
	
	k = 0;		
	for(i = 0; i < npars; i++){
	  for(j = 0; j < nitems; j++){
	    Preturn[k] = 0;
	    k++;
		}
	}		
	
	//BIG LOOP OVER ITEMS 		
	for(int item = 0; item < nitems; item++)
	{
		//Load the values		
		for(i = 0; i < npars; i++){
	   StepLimit[i] = 0.5;
	   LastL[i] = 0;
	  }	
	  for(i = 0; i < nquad; i++){
		  r1[i] = fullr1[item][i];
		  N[i] = fullN[item][i];
		}
		for(i = 0; i < nfact; i++)
		  a[i] = fullpars[item][i];
	  d = fullpars[item][nfact];
		g = guess[item]; 	  
	  for(int loop = 0; loop < 100; loop++)//MAX LOOP	  			  		    
	  {	
		  itemtrace(P, a, &d, PTheta, &g, &nfact, &nquad);	  	  		  
		  for(i = 0; i < nquad; i++){      
		  	PQ[i] = P[i] * (1 - P[i]) * Pprior[i];
		  	DIF[i] = ((r1[i] / N[i]) - P[i]) * Pprior[i];	  	
		  }		  
		  //gradient		  
		  for(i = 0 ; i < nfact; i++){
		  	for(j = 0; j < nquad; j++)
		      tempTheta[j] = Theta[j][i]; 
		    arrayprod3(temparray, N, DIF, tempTheta, &nquad);
		    L[i] = arraysum(temparray, &nquad);		     
		  }		  	    
		  arrayprod2(temparray, N, DIF, &nquad);
		  L[nfact] = arraysum(temparray, &nquad);		       
		  //hessian
		  for(i = 0 ; i < nfact; i++){
		    for(j = 0; j < nquad; j++)
		      tempTheta[j] = Theta[j][i] * Theta[j][i];
		    arrayprod3(temparray, N, PQ, tempTheta, &nquad);
		    LL[i][i] = (-1.0) * arraysum(temparray, &nquad);
		  }	  	    
		  arrayprod2(temparray, N, PQ, &nquad);	  
		  LL[nfact][nfact] = (-1.0) * arraysum(temparray, &nquad);		  
		  for(i = 0; i < nfact; i++){
		  	for(j = 0; j < nfact; j++){		      		  	  		  
		  	  if( i < j ){
		  	  	for(k = 0; k < nquad; k++)	  	    
		          tempTheta[k] = Theta[k][i] * Theta[k][j];
		  	  	arrayprod3(temparray, N, PQ, tempTheta, &nquad);
		  	  	LL[i][j] = LL[j][i] = (-1.0) * arraysum(temparray, &nquad);
		  	  }
		  	}
		  }
		  for(i = 0 ; i < nfact; i++){
		  	for(j = 0; j < nquad; j++)
		      tempTheta[j] = Theta[j][i];
		    arrayprod3(temparray, N, PQ, tempTheta, &nquad);
		    LL[nfact][i] = LL[i][nfact] = (-1.0) * arraysum(temparray, &nquad);
		  }		  		  
		  //priors
		  if(fullparprior[item][0] > 1.0){	  
			  d2 = 1;
		  	for(i = 0; i < nfact; i++)
		  	  d2 += a[i] * a[i];	  	  
		  	c = 2.0*(fullparprior[item][0] - 1) / d2*d2;
		  	for(i = 0; i < nfact; i++)
		  	  L[i] -= c * a[i];
		  	for(i = 0; i < nfact; i++)
		  	  betaprior[i][i] = d2 - 2*a[i]*a[i];
		  	for(i = 0; i < nfact; i++)
		  	  for(j = 0; j < nfact; j++)
		  	    if(i < j) betaprior[i][j] = betaprior[j][i] = -2*a[i]*a[j];
		  	for(i = 0; i < nfact; i++)
		  	  for(j = 0; j < nfact; j++)
		  	    LL[i][j] += betaprior[i][j];
	  	}        	 	  		    
		  if(fullparprior[item][2] > 0.0){
		  	normprior = dnorm(fullparprior[item][1],fullparprior[item][1], fullparprior[item][2],0) -
		  	  dnorm(d, fullparprior[item][1], fullparprior[item][2], 0);		  
			  L[nfact] -=  2*normprior;			    		  			  			  	
			  for(i = 0; i < npars; i++){
			    LL[i][nfact] -= 2*normprior;
			    LL[nfact][i] = LL[i][nfact];
			  }
		  }
	    //Invert and NR correction
	    for(i = 0; i < npars; i++) 		  
		    for(j = 0; j < npars; j++)	  
		      LLinv[i][j] = LL[i][j];  		  
		  F77_CALL(dgetrf)(&npars, &npars, &LLinv[0][0], &npars, ipiv, &info);
		  if(info > 0){
		  	for(i = 0; i < npars; i++)
		  	  LLinv[i][i] += .0001;
		    F77_CALL(dgetrf)(&npars, &npars, &LLinv[0][0], &npars, ipiv, &info);	
		  }	
		  F77_CALL(dgetri)(&npars, &LLinv[0][0], &npars, ipiv, work, &npars, &info); 		  		  		  		  
		  for(i = 0; i < npars; i++)
		    correction[i] = 0;
		  for(i = 0; i < npars; i++)
		    for(j = 0; j < npars; j++)
		      correction[i] += LLinv[i][j] * L[j];		      
		  //stop condition		   		  
		  for(i = 0; i < npars; i++){
		  	corSign[i] = 1.0;
		  	if(correction[i] < 0) corSign[i] = -1.0;
		  }		  
		  temp = 0.0; 
		  for(i = 0; i < npars; i++){		  	
		    temp += corSign[i]*correction[i];
		  }       	  
		  if(temp < .00001) break;		  	
		  //rate checking	  		  	
		  for(i = 0; i < npars; i++){
		    if((corSign[i] * correction[i]) > StepLimit[i])		   	  
		      correction[i] = corSign[i] * StepLimit[i];
		  }    		  
		  for(i = 0; i < nfact; i++) 
		    a[i] -= correction[i];
		  d -= correction[nfact];     		
		  for(i = 0; i < nfact; i++){		    
		    if(L[i] * LastL[i] < 0.0){
		    	a[i] += 0.5 * correction[i];
		    	StepLimit[i] = 0.5 * StepLimit[i];
		    }
		  }
		  if(L[nfact] * LastL[nfact] < 0.0){
		    d += 0.5 * correction[nfact];
		    StepLimit[nfact] = 0.5 * StepLimit[nfact];
		  }		  
		  for(i = 0; i < npars; i++)
		    LastL[i] = L[i];
		  iter++;		   
		} //END MAX LOOP		
		for(i = 0; i < nfact; i++)
		  fullpars[item][i] = a[i];
	    fullpars[item][nfact] = d;
			 
  }//  END BIG LOOP						
	k = 0;		
	for(i = 0; i < npars; i++){
	  for(j = 0; j < nitems; j++){
	    Preturn[k] = fullpars[j][i];
	    k++;
		}
	}		
	UNPROTECT(8);		
	return(Rreturn);	
}
