#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>

static void polyOuter(double *d2Louter, const double *PThetas, const double *Pk,
	const double *Pk_1, const double *PQ_1,	const double *PQ, 
	const double *dif1sq, const double *dif1, const unsigned int *nfact, 
	const unsigned int *N){
	
	unsigned int i, j, n;
	double out[*nfact][*nfact], Thetas[*N][*nfact], outer[*nfact][*nfact],
		temp[*nfact];

	for(i = 0; i < *nfact; i++)
		for(j = 0; j < *nfact; j++)
			out[i][j] = 0.0;
	for(i = 0; i < *nfact; i++)
		for(n = 0; n < *N; n++)
			Thetas[n][i] = PThetas[n + (*N)*i];	
	for(n = 0; n < *N; n++){
		for(i = 0; i < *nfact; i++)
			for(j = 0; j < *nfact; j++)
				outer[i][j] = Thetas[n][i] * Thetas[n][j];
		for(i = 0; i < *nfact; i++)
			temp[i] =  (PQ_1[n] * Thetas[n][i] - PQ[n] * Thetas[n][i]);
		for(i = 0; i < *nfact; i++)
			for(j = 0; j < *nfact; j++)				
				out[i][j] += (-1) * dif1sq[n] * temp[i] * temp[j] +  
				(dif1[n] * (Pk_1[n] * (1.0 - Pk_1[n]) * (1.0 - 2.0 * Pk_1[n]) * 
				outer[i][j] - Pk[n] * (1.0 - Pk[n]) * (1.0 - 2.0 * Pk[n]) * outer[i][j]));
	}				
	n = 0;
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *nfact; j++){
			d2Louter[n] = out[j][i];	
			n++;
		}
	}	
}

static void itemtrace(double *P, const double *a, 
  const double *d, const double *PTheta, const double *g, 
  const unsigned int *nfact, const unsigned int *nquad)
{	
	double z[*nquad];
	int i, j, loc[*nfact];

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

static void Prob(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g)
{
	double Ps[*N][*k + 1], Pdif[*N][*k], p1[*N], tmp;
	unsigned int i, m = 0;
	int j;

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
	for(j = 0; j < *k; j++){
		for(i = 0; i < *N; i++){
			if(Pdif[i][j] < .00000001) Pdif[i][j] = .00000001;
			if(*k == 2) Pdif[i][j] = 1 - Pdif[i][j];
			P[m] = Pdif[i][j];
			m++;
		}
	}
}

static void ProbComp(double *P, const unsigned int *k, const unsigned int *N, 
	const unsigned int *nfact, const double *theta, const double *a, 
	const double *d, const double *g)
{
	double Theta[*N], tmp[*N], zerog = 0.0, tmpa, tmpd;
	unsigned int i, j, onenfact = 1;	
	for(j = 0; j < *N; j++)
		tmp[j] = 1.0;
	for(i = 0; i < *nfact; i++){
		for(j = 0; j < *N; j++)
			Theta[j] = theta[j + i*(*N)];
		tmpa = a[i];
		tmpd = d[i];
		itemtrace(P, &tmpa, &tmpd, Theta, &zerog, &onenfact, N);
		for(j = 0; j < *N; j++)
			tmp[j] *= P[j];
	}
	for(j = 0; j < *N; j++){
		P[j] = *g + (1 - *g)*tmp[j];
		P[j + *N] = 1.0 - P[j];
	}
}


