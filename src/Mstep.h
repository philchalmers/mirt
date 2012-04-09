#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

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


