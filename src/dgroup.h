#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

static void matrixMult(double *c, const double *a, const double *b, 
	const unsigned int *dim)
{
	double A[*dim][*dim], B[*dim][*dim], C[*dim][*dim];
	unsigned int i, j, k = 0;

	for (j = 0; j < *dim; j++){ 
		for (i = 0; i < *dim; i++){ 		
			A[i][j] = a[k];
			k++;
		}
	}
	k = 0;
	for (j = 0; j < *dim; j++){ 
		for (i = 0; i < *dim; i++){ 		
			B[i][j] = b[k];
			k++;
		}
	}
	for (i = 0; i < *dim; i++){ 
		for (j = 0; j < *dim; j++) {
			C[i][j] = 0;
			for (k = 0; k < *dim; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
	k = 0;
	for (j = 0; j < *dim; j++) {
		for (i = 0; i < *dim; i++){ 		
			c[k] = C[i][j]; 
			k++;
		}
	}   
}

static void matrixMult4(double *e, const double *a, const double *b,
	const double *c, const double *d, const unsigned int *dim)
{
	double tmp1[*dim * (*dim)], tmp2[*dim * (*dim)];
	matrixMult(tmp1, a, b, dim);
	matrixMult(tmp2, tmp1, c, dim);
	matrixMult(e, tmp2, d, dim);
}


static double tr(double *a, const unsigned int *dim)
{	
	double trace = 0.0;
	unsigned int i, j, k = 0;

	for(j = 0; j < *dim; j++){
		for(i = 0; i < *dim; i++){
			if(i == j)
				trace += a[k];			
			k++;
		}
	}	
	return trace;
}

static void matrixSub(double *c, const double *a,const  double *b, 
	const unsigned int *dim)
{	
	unsigned int i;
	for(i = 0; i < *dim*(*dim); i++)		
		c[i] = a[i] - b[i];
}

static void outer(double *c, const double *a, const double *b, 
	const  unsigned int *dim)
{
	unsigned int i, j, k = 0;
	for(i = 0; i < *dim; i++){
		for(j = 0; j < *dim; j++){
			c[k] = a[j] * b[i];
			k++;
		}
	}
}

static double inner(double *a, const double *b, const double *c, 
	const unsigned int *dim)
{
	unsigned int i, j, k = 0;
	double tmp[*dim], B[*dim][*dim], ret = 0.0;
				
	for(i = 0; i < *dim; i++){
		tmp[i] = 0.0;
		for(j = 0; j < *dim; j++){
			B[j][i] = b[k];
			k++;
		}
	}		
	for(i = 0; i < *dim; i++){
		for(j = 0; j < *dim; j++){
			tmp[i] += a[j] * B[j][i];
			k++;
		}
	}
	for(i = 0; i < *dim; i++)
		ret += tmp[i] * c[i];
	return ret;
}

static void symMat(double *dsig, const unsigned int *nfact)
{
	unsigned int i, j, k = 0;
	double tmp[*nfact][*nfact];
	
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *nfact; j++){
			tmp[i][j] = dsig[k];
			k++;
		}
	}
	for(i = 0; i < *nfact; i++)
		for(j = 0; j < *nfact; j++)
			if(i < j)
				tmp[j][i] = tmp[i][j];
	k = 0;
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *nfact; j++){
			dsig[k] = tmp[i][j];
			k++;
		}
	}	
}

