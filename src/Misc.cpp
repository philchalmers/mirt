#include"Misc.h"

NumericMatrix polyOuter(NumericMatrix Thetas, NumericVector Pk,
	NumericVector Pk_1, NumericVector PQ_1, NumericVector PQ, 
	NumericVector dif1sq, NumericVector dif1)
{
	int i, j, n, nfact, N;
	nfact = Thetas.ncol();
	N = Thetas.nrow();
	NumericMatrix d2Louter(nfact,nfact), outer(nfact,nfact);
	NumericVector temp(nfact);
	d2Louter.fill(0.0);
	
	for(n = 0; n < N; n++){
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				outer(i,j) = Thetas(n,i) * Thetas(n,j);
		for(i = 0; i < nfact; i++)
			temp(i) =  (PQ_1(n) * Thetas(n,i) - PQ(n) * Thetas(n,i));
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)				
				d2Louter(i,j) += (-1) * dif1sq(n) * temp(i) * temp(j) +  
				    (dif1(n) * (Pk_1(n) * (1.0 - Pk_1(n)) * (1.0 - 2.0 * Pk_1(n)) * 
				    outer(i,j) - Pk(n) * (1.0 - Pk(n)) * (1.0 - 2.0 * Pk(n)) * outer(i,j)));
	}
	return d2Louter;		
}

NumericVector itemTrace(NumericVector a, const double *d, 
        NumericMatrix Theta, const double *g)
{	
	int i, j;
    int nquad = Theta.nrow();
    int nfact = Theta.ncol();
	NumericVector P(nquad), z(nquad);

	for (i = 0; i <	nquad; i++){
	    z(i) = 0.0;
		for (j = 0; j <	nfact; j++){		
			z(i) += 1.702 * a(j) * Theta(i,j);  		
		}
		z(i) += *d * 1.702;
	}	
	for (i = 0; i < nquad; i++) 
		P(i) = *g + (1 - *g) * (exp(z(i))/(1 + exp(z(i))));
	
	return P;		
}

NumericMatrix Prob(NumericMatrix Theta, NumericVector a,
        NumericVector zetas, const double *g)
{
	int i, j;
    int N = Theta.nrow();
    int k = zetas.length() + 1;
    double tmp;
	NumericVector p1(N);
	NumericMatrix Ps(N,k+1), Pdif(N,k), P(N,k);

	for(i = 0; i < N; i++){
		Ps(i,0) = 1.0;
		Ps(i,k) = 0.0;
	}
	for(j = 0; j < (k - 1); j++){
		tmp = zetas(j);
		p1 = itemTrace(a, &tmp, Theta, g);
		for(i = 0; i < N; i++)
			Ps(i,j+1) = p1(i);
	}
	for(j = (k - 1); j >= 0; j--)
		for(i = 0; i < N; i++)
			Pdif(i,j) = Ps(i,j) - Ps(i,j+1);				
	for(j = 0; j < k; j++){
		for(i = 0; i < N; i++){
			if(Pdif(i,j) < .00000001) Pdif(i,j) = .00000001;
			if(k == 2) Pdif(i,j) = 1.0 - Pdif(i,j);
			P(i,j) = Pdif(i,j);
		}
	}
	return P;
}

NumericMatrix ProbComp(NumericMatrix Theta, NumericVector a, 
        NumericVector zetas, const double *g)
{
	int i, j;
    int nfact = Theta.ncol();
    int N = Theta.nrow();
	NumericMatrix Pret(N,2);    
	NumericVector P(N), p1(N), tmpa(1);
	double zerog = 0.0, tmpd;
	P.fill(1.0);

	for(j = 0; j < nfact; j++){
		tmpa(0) = a(j);
		tmpd = zetas(j);
		p1 = itemTrace(tmpa, &tmpd, Theta, &zerog);
		for(i = 0; i < N; i++)
			P(i) *= p1(i);
	}
	for(i = 0; i < N; i++){
		Pret(i,0) = *g + (1.0 - *g) * P(i);
		Pret(i,1) = 1.0 - Pret(i,0);
	}
	return Pret;
}


