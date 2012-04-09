#include"Misc.h"

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




