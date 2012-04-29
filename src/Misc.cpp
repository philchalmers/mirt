#include"Misc.h"

void polyOuter(double *d2Louter, NumericMatrix Thetas, const double *Pk,
	const double *Pk_1, const double *PQ_1,	const double *PQ, 
	const double *dif1sq, const double *dif1)
{
	int i, j, n, nfact, N;
	nfact = Thetas.ncol();
	N = Thetas.nrow();
	double out[nfact][nfact], outer[nfact][nfact], temp[nfact];

	for(i = 0; i < nfact; i++)
		for(j = 0; j < nfact; j++)
			out[i][j] = 0.0;
	for(n = 0; n < N; n++){
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				outer[i][j] = Thetas(n,i) * Thetas(n,j);
		for(i = 0; i < nfact; i++)
			temp[i] =  (PQ_1[n] * Thetas(n,i) - PQ[n] * Thetas(n,i));
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)				
				out[i][j] += (-1) * dif1sq[n] * temp[i] * temp[j] +  
				    (dif1[n] * (Pk_1[n] * (1.0 - Pk_1[n]) * (1.0 - 2.0 * Pk_1[n]) * 
				    outer[i][j] - Pk[n] * (1.0 - Pk[n]) * (1.0 - 2.0 * Pk[n]) * outer[i][j]));
	}				
	n = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nfact; j++){
			d2Louter[n] = out[j][i];	
			n++;
		}
	}	
}

void itemtrace(double *P, const double *a, const double *d, 
        NumericMatrix Theta, const double *g)
{	
	int i, j;
    int nquad = Theta.nrow();
    int nfact = Theta.ncol();
	double z[nquad];

	for (i = 0; i <	nquad; i++){
	    z[i] = 0.0;
		for (j = 0; j <	nfact; j++){		
			z[i] += 1.702 * a[j] * Theta(i,j);  		
		}
		z[i] += *d * 1.702;
	}	
	for (i = 0; i < nquad; i++) 
		P[i] = *g + (1 - *g) * (exp(z[i])/(1 + exp(z[i])));		
}

NumericMatrix Prob(NumericMatrix Theta, const double *a,
        NumericVector zetas, const double *g)
{
	int i, j;
    int N = Theta.nrow();
    int k = zetas.length() + 1;
    double Ps[N][k+1], Pdif[N][k], p1[N], tmp;
	NumericMatrix P(N,k);

	for(i = 0; i < N; i++){
		Ps[i][0] = 1.0;
		Ps[i][k] = 0.0;
	}
	for(j = 0; j < (k - 1); j++){
		tmp = zetas[j];
		itemtrace(p1, a, &tmp, Theta, g);
		for(i = 0; i < N; i++)
			Ps[i][j+1] = p1[i];
	}
	for(j = (k - 1); j >= 0; j--)
		for(i = 0; i < N; i++)
			Pdif[i][j] = Ps[i][j] - Ps[i][j+1];				
	for(j = 0; j < k; j++){
		for(i = 0; i < N; i++){
			if(Pdif[i][j] < .00000001) Pdif[i][j] = .00000001;
			if(k == 2) Pdif[i][j] = 1.0 - Pdif[i][j];
			P(i,j) = Pdif[i][j];
		}
	}
	return P;
}

NumericMatrix ProbComp(NumericMatrix Theta, const double *a, 
        NumericVector zetas, const double *g)
{
	int i, j;
    int nfact = Theta.ncol();
    int N = Theta.nrow();
	NumericMatrix P(N,1), Pret(N,2);
    P.fill(1.0);
	double p1[N], zerog = 0.0, tmpa, tmpd;

	for(j = 0; j < nfact; j++){
		tmpa = a[j];
		tmpd = zetas[j];
		itemtrace(p1, &tmpa, &tmpd, Theta, &zerog);
		for(i = 0; i < N; i++)
			P[i] *= p1[i];
	}
	for(i = 0; i < N; i++){
		Pret(i,0) = *g + (1.0 - *g) * P[i];
		Pret(i,1) = 1.0 - Pret(i,0);
	}
	return Pret;
}


