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

SEXP dparsPoly(SEXP Rprob, SEXP RThetas, SEXP Rdat, SEXP Rnzeta, 
	SEXP Rnfact, SEXP RN){
		
	unsigned int i, j, k, nzeta, nfact, N; 
	int *Pdat;
	double *Pprob, *PThetas;		
	
	PROTECT(Rprob = AS_NUMERIC(Rprob));	
	PROTECT(RThetas = AS_NUMERIC(RThetas));
	PROTECT(Rdat = AS_INTEGER(Rdat));
	PROTECT(Rnzeta = AS_INTEGER(Rnzeta));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RN = AS_INTEGER(RN));
	Pprob = NUMERIC_POINTER(Rprob);
	PThetas = NUMERIC_POINTER(RThetas);
	Pdat = INTEGER_POINTER(Rdat);
	nzeta = NUMERIC_VALUE(Rnzeta);
	nfact = NUMERIC_VALUE(Rnfact);
	N = NUMERIC_VALUE(RN);	

	double P[N][nzeta+2], PQfull[N][nzeta+2], Thetas[N][nfact], dL[nfact + nzeta],
		d2L[nfact + nzeta][nfact + nzeta], Pk[N], Pk_1[N], Pk_p1[N], PQ_1[N], PQ[N], 
		PQ_p1[N], Pk_1Pk[N], Pk_Pkp1[N], dif1[N], dif1sq[N], dif2[N], dif2sq[N],
		tmp1[N], tmp2[N], tmp3[N], csums[nfact], mattmp[N][nfact], d2Louter[nfact * nfact],
		tmp;
	unsigned int factind[nfact], dat[N][nzeta + 1];
	k = 0;
	for(j = 0; j < (nzeta + 2); j++){
		for(i = 0; i < N; i++){
			P[i][j] = Pprob[k];
			PQfull[i][j] = Pprob[k] * (1.0 - Pprob[k]);
			k++;
		}
	}
	k = 0;
	for(j = 0; j < (nzeta + 1); j++){
		for(i = 0; i < N; i++){			
			dat[i][j] = Pdat[k];
			k++;
		}
	}
	k = 0;
	for(j = 0; j < nfact; j++){
		factind[j] = nzeta + j;
		for(i = 0; i < N; i++){
			Thetas[i][j] = PThetas[k];
			k++;
		}
	}
	for(j = 0; j < (nfact + nzeta); j++){	
		dL[j] = 0.0;
		for(i = 0; i < (nfact + nzeta); i++)
			d2L[i][j] = 0.0;
	}		

	for(j = 0; j < (nzeta + 1); j++){
		if(j < nzeta){
			for(i = 0; i < N; i++){
				Pk_1[i] = P[i][j];
				Pk[i] = P[i][j + 1];
				Pk_p1[i] = P[i][j + 2];
				PQ_1[i] = PQfull[i][j];
				PQ[i] = PQfull[i][j + 1];
				PQ_p1[i] = PQfull[i][j + 2];
				Pk_1Pk[i] = Pk_1[i] - Pk[i];
				Pk_Pkp1[i] = Pk[i] - Pk_p1[i];
				if(Pk_1Pk[i] < 1e-10) Pk_1Pk[i] = 1e-10;
				if(Pk_Pkp1[i] < 1e-10) Pk_Pkp1[i] = 1e-10;
				dif1[i] = dat[i][j] / Pk_1Pk[i];
				dif1sq[i] = dat[i][j] / (Pk_1Pk[i] * Pk_1Pk[i]);
				dif2[i] = dat[i][j+1] / Pk_Pkp1[i];
				dif2sq[i] = dat[i][j+1] / (Pk_Pkp1[i] * Pk_Pkp1[i]);
			}			
			tmp = 0.0;
			for(i = 0; i < N; i++)
				tmp += (-1.0) * PQ[i] * (dif1[i] - dif2[i]);			
			dL[j] = tmp;			
			tmp = 0.0;
			for(i = 0; i < N; i++)
				tmp += (-1.0) * PQ[i] * PQ[i] * (dif1sq[i] + dif2sq[i]) -				
					(dif1[i] - dif2[i]) * (Pk[i] * (1.0 - Pk[i]) * (1.0 - 2.0*Pk[i]));			
			d2L[j][j] = tmp;
			if(j < (nzeta - 1)){
				tmp = 0.0;
				for(i = 0; i < N; i++)
					tmp += dif2sq[i] * PQ_p1[i] * PQ[i];
				d2L[j][j + 1] = tmp;
				d2L[j+1][j] = tmp;
			}
			for(i = 0; i < N; i++){
				tmp1[i] = (-1.0) * dif2sq[i] * PQ[i] * (PQ[i] - PQ_p1[i]);
				tmp2[i] = dif1sq[i] * PQ[i] * (PQ_1[i] - PQ[i]);
				tmp3[i] = (dif1[i] - dif2[i]) * (Pk[i] * (1.0 - Pk[i]) * (1.0 - 2.0*Pk[i]));
			}
			for(k = 0; k < nfact; k++){
				csums[k] = 0.0;
				for(i = 0; i < N; i++){
					mattmp[i][k] = tmp1[i] * Thetas[i][k] + tmp2[i] * Thetas[i][k] - 
						tmp3[i] * Thetas[i][k];
					csums[k] += mattmp[i][k];
				}
			}
			for(i = 0; i < nfact; i++){
				d2L[j][factind[i]] = csums[i];
				d2L[factind[i]][j] = csums[i];
			}			
		} else {					
			for(i = 0; i < N; i++){
				Pk_1[i] = P[i][j];
				Pk[i] = P[i][j + 1];			
				PQ_1[i] = PQfull[i][j];
				PQ[i] = PQfull[i][j + 1];			
				Pk_1Pk[i] = Pk_1[i] - Pk[i];			
				if(Pk_1Pk[i] < 1e-10) Pk_1Pk[i] = 1e-10;			
				dif1[i] = dat[i][j] / Pk_1Pk[i];
				dif1sq[i] = dat[i][j] / (Pk_1Pk[i] * Pk_1Pk[i]);			
			}	
		}
		for(k = 0; k < nfact; k++){
			csums[k] = 0.0;
			for(i = 0; i < N; i++){
				mattmp[i][k] = dif1[i] * (PQ_1[i] - PQ[i]) * Thetas[i][k];
				csums[k] += mattmp[i][k];
			}
		}
		for(i = 0; i < nfact; i++)
			dL[factind[i]] += csums[i];			
		
		polyOuter(d2Louter, PThetas, Pk, Pk_1, PQ_1, PQ, dif1sq, dif1, &nfact, &N);
		unsigned int m=0;
		for(k = 0; k < nfact; k++){
			for(i = 0; i < nfact; i++){
				d2L[factind[i]][factind[k]] += d2Louter[m];
				m++;
			}
		}		
	}

	SEXP Rgrad, Rhess, list_names, list;
	PROTECT(Rgrad = allocVector(REALSXP,nfact + nzeta));	
	PROTECT(Rhess = allocMatrix(REALSXP,nfact + nzeta,nfact + nzeta));	
	for (i = 0; i < (nfact + nzeta); i++)
		REAL(Rgrad)[i] = dL[i];
	k = 0;
	for (j = 0; j < (nfact + nzeta); j++){
		for (i = 0; i < (nfact + nzeta); i++){
			REAL(Rhess)[k] = d2L[i][j];
			k++;
		}
	}	
		
	//set list names		
	char *names[2] = {"grad","hess"};
	PROTECT(list_names = allocVector(STRSXP,2));
	for(i = 0; i < 2; i++)
		SET_STRING_ELT(list_names, i, mkChar(names[i]));
	//set list
	PROTECT(list = allocVector(VECSXP,2));
	SET_VECTOR_ELT(list, 0, Rgrad);
	SET_VECTOR_ELT(list, 1, Rhess);	
	setAttrib(list, R_NamesSymbol, list_names); 
		
	UNPROTECT(10);		
	return(list);
}

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


SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP Rnfact, SEXP RN){
	
	SEXP Rreturn;			
	unsigned int i, j, n, nfact, N;
	double *Preturn,*PThetas,*PQ;	
	
	//set pointers and protect
	PROTECT(RThetas = AS_NUMERIC(RThetas));	
	PROTECT(RPQ = AS_NUMERIC(RPQ));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(RN = AS_INTEGER(RN));
	PThetas = NUMERIC_POINTER(RThetas);
	PQ = NUMERIC_POINTER(RPQ);
	nfact = NUMERIC_VALUE(Rnfact);
	N = NUMERIC_VALUE(RN);
	PROTECT(Rreturn = allocMatrix(REALSXP,nfact,nfact));		
	Preturn = NUMERIC_POINTER(Rreturn);

	double out[nfact][nfact], Thetas[N][nfact];
	for(i = 0; i < nfact; i++)
		for(j = 0; j < nfact; j++)
			out[i][j] = 0.0;
	for(i = 0; i < nfact; i++)
		for(n = 0; n < N; n++)
			Thetas[n][i] = PThetas[n + N*i];	
	for(n = 0; n < N; n++)
		for(i = 0; i < nfact; i++)
			for(j = 0; j < nfact; j++)
				out[i][j] += Thetas[n][i] * Thetas[n][j] * PQ[n];
	n = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < nfact; j++){
			Preturn[n] = out[j][i];	
			n++;
		}
	}
		
	UNPROTECT(5);	
	return(Rreturn);
}

SEXP drawThetas(SEXP Runif, SEXP Rden0, SEXP Rden1, SEXP Rlambdas, SEXP Rzetas, 
	SEXP Rguess, SEXP Rtheta0, SEXP Rtheta1, SEXP Rfulldata, SEXP Ritemloc,
	SEXP RK, SEXP RJ, SEXP RN, SEXP Rnfact){

	SEXP Rreturn;
	unsigned int i, j, k, m, nfact, J, N, Ksums = 0, max = 2;
	int *itemloc,*K,*Pfulldata;
	double *Preturn,*Plambdas,*zetas,*guess,*Ptheta0,*Ptheta1,*unif, 
		*den0,*den1;
		
	PROTECT(Runif = AS_NUMERIC(Runif));
	PROTECT(Rden0 = AS_NUMERIC(Rden0));
	PROTECT(Rden1 = AS_NUMERIC(Rden1));
	PROTECT(Rlambdas = AS_NUMERIC(Rlambdas));
	PROTECT(Rzetas = AS_NUMERIC(Rzetas));
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(Rtheta0 = AS_NUMERIC(Rtheta0));
	PROTECT(Rtheta1 = AS_NUMERIC(Rtheta1));
	PROTECT(Rfulldata = AS_INTEGER(Rfulldata));
	PROTECT(Ritemloc = AS_INTEGER(Ritemloc));
	PROTECT(RK = AS_INTEGER(RK));	
	PROTECT(RJ = AS_INTEGER(RJ));	
	PROTECT(RN = AS_INTEGER(RN));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	unif = NUMERIC_POINTER(Runif);
	den0 = NUMERIC_POINTER(Rden0);
	den1 = NUMERIC_POINTER(Rden1);
	Plambdas = NUMERIC_POINTER(Rlambdas);
	zetas = NUMERIC_POINTER(Rzetas);
	guess = NUMERIC_POINTER(Rguess);
	Ptheta0 = NUMERIC_POINTER(Rtheta0);
	Ptheta1 = NUMERIC_POINTER(Rtheta1);
	Pfulldata = INTEGER_POINTER(Rfulldata);
	itemloc = INTEGER_POINTER(Ritemloc);	
	K = INTEGER_POINTER(RK);	
	J = NUMERIC_VALUE(RJ);
	N = NUMERIC_VALUE(RN);
	nfact = NUMERIC_VALUE(Rnfact);	
	for(i = 0; i < J; i++){
		Ksums += K[i]; 
		if(max < K[i]) max = K[i];
	}
	
	PROTECT(Rreturn = NEW_NUMERIC(N + 1));
	Preturn = NUMERIC_POINTER(Rreturn);
	double a[nfact], d[max], g, lambdas[J][nfact], irt0[N], irt1[N], 
		accept[N], Plong_1[N * max], Plong_0[N * max], cdloglik = 0;
	unsigned int loc = 0, location[J];
	
	k = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < J; j++){
			lambdas[j][i] = Plambdas[k];
			k++;
		}
	}		
	for(i = 0; i < J; i++)		
		location[i] = itemloc[i] * N;	
	for(i = 0; i < N; i++){
		irt0[i] = 0.0;
		irt1[i] = 0.0;
	}	
	for(unsigned int item = 0; item < J; item++){		
		k = K[item];
		for(i = 0; i < nfact; i++)
			a[i] = lambdas[item][i];		
		for(i = 0; i < (k-1); i++) 
			d[i] = zetas[i + loc];
		g = guess[item];		
		Prob(Plong_0, &k, &N, &nfact, Ptheta0, a, d, &g);			
		Prob(Plong_1, &k, &N, &nfact, Ptheta1, a, d, &g);		
		m = 0;
		for(j = 0; j < k; j++){
			for(i = 0; i < N; i++){				
				if(Pfulldata[m + location[item]]){
					irt0[i] += log(Plong_0[m]);
					irt1[i] += log(Plong_1[m]);
				}				
				m++;
			}
		}		
		loc += k - 1;
	}	
	for(i = 0; i < N; i++){		
		irt0[i] += log(den0[i]);
		irt1[i] += log(den1[i]);
		accept[i] = irt1[i] - irt0[i];		
		if(accept[i] > 0.0) accept[i] = 0.0;
		if(unif[i] < exp(accept[i])) accept[i] = 1.0;
			else accept[i] = 0.0;
		Preturn[i] = accept[i];
	}	
	for(i = 0; i < N; i++){		
		if(accept[i]) cdloglik += irt1[i];
		else cdloglik += irt0[i];
	}
	Preturn[N] = cdloglik;
	
	UNPROTECT(15);	
	return(Rreturn);	
}

SEXP logLik(SEXP Rlambdas, SEXP Rzetas, SEXP Rguess, SEXP Rtheta0,
	SEXP Rfulldata, SEXP Ritemloc, SEXP RK, SEXP RJ, SEXP RN, SEXP Rnfact){

	SEXP Rreturn;
	unsigned int i, j, k, m, nfact, J, N, Ksums = 0, max = 2;
	int *itemloc,*K,*Pfulldata;
	double *Preturn,*Plambdas,*zetas,*guess,*Ptheta0;
		
	PROTECT(Rlambdas = AS_NUMERIC(Rlambdas));
	PROTECT(Rzetas = AS_NUMERIC(Rzetas));
	PROTECT(Rguess = AS_NUMERIC(Rguess));
	PROTECT(Rtheta0 = AS_NUMERIC(Rtheta0));	
	PROTECT(Rfulldata = AS_INTEGER(Rfulldata));
	PROTECT(Ritemloc = AS_INTEGER(Ritemloc));
	PROTECT(RK = AS_INTEGER(RK));	
	PROTECT(RJ = AS_INTEGER(RJ));	
	PROTECT(RN = AS_INTEGER(RN));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));		
	Plambdas = NUMERIC_POINTER(Rlambdas);
	zetas = NUMERIC_POINTER(Rzetas);
	guess = NUMERIC_POINTER(Rguess);
	Ptheta0 = NUMERIC_POINTER(Rtheta0);	
	Pfulldata = INTEGER_POINTER(Rfulldata);
	itemloc = INTEGER_POINTER(Ritemloc);	
	K = INTEGER_POINTER(RK);	
	J = NUMERIC_VALUE(RJ);
	N = NUMERIC_VALUE(RN);
	nfact = NUMERIC_VALUE(Rnfact);	
	for(i = 0; i < J; i++){
		Ksums += K[i]; 
		if(max < K[i]) max = K[i];
	}
	
	PROTECT(Rreturn = NEW_NUMERIC(N));
	Preturn = NUMERIC_POINTER(Rreturn);
	double a[nfact], d[max], g, lambdas[J][nfact], irt0[N],
		Plong_0[N * max], cdloglik;
	unsigned int loc = 0, location[J];
	
	k = 0;
	for(i = 0; i < nfact; i++){
		for(j = 0; j < J; j++){
			lambdas[j][i] = Plambdas[k];
			k++;
		}
	}		
	for(i = 0; i < J; i++)		
		location[i] = itemloc[i] * N;	
	for(i = 0; i < N; i++)
		irt0[i] = 1.0;		
		
	for(unsigned int item = 0; item < J; item++){		
		k = K[item];
		for(i = 0; i < nfact; i++)
			a[i] = lambdas[item][i];		
		for(i = 0; i < (k-1); i++) 
			d[i] = zetas[i + loc];
		g = guess[item];		
		Prob(Plong_0, &k, &N, &nfact, Ptheta0, a, d, &g);							
		m = 0;
		for(j = 0; j < k; j++){
			for(i = 0; i < N; i++){				
				if(Pfulldata[m + location[item]])
					irt0[i] *= Plong_0[m];													
				m++;
			}
		}		
		loc += k - 1;
	}	
	for(i = 0; i < N; i++){				 
		cdloglik = irt0[i];
		if(cdloglik < .000000000001) cdloglik = .000000000001;
		Preturn[i] = cdloglik;
	}
	
	UNPROTECT(11);	
	return(Rreturn);	
}
