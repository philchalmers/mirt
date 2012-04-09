#include"dgroup.h"

SEXP dgroup(SEXP RinvSig, SEXP RcMeans,	SEXP RZdif, SEXP RN, SEXP Rnfact, 
        SEXP Rnpars) 
{   
	//SEXP Rreturn;			
	unsigned int i, j, k, N, nfact, npars, nsig;	
	double *invSig, *cMeans, *Zdif;
	
	PROTECT(RinvSig = AS_NUMERIC(RinvSig));
	PROTECT(RcMeans = AS_NUMERIC(RcMeans));
	PROTECT(RZdif = AS_NUMERIC(RZdif));
	PROTECT(RN = AS_INTEGER(RN));
	PROTECT(Rnfact = AS_INTEGER(Rnfact));
	PROTECT(Rnpars = AS_INTEGER(Rnpars));
	invSig = NUMERIC_POINTER(RinvSig);
	cMeans = NUMERIC_POINTER(RcMeans);
	Zdif = NUMERIC_POINTER(RZdif);
	N = INTEGER_VALUE(RN);
	nfact = INTEGER_VALUE(Rnfact);
	npars = INTEGER_VALUE(Rnpars);
	nsig = npars - nfact;
		 
	double derv1[npars], derv2[npars], du1[nfact], du2[nfact], dsig1[nsig],
		dsig2[nsig], dZ[nsig], dinvSig2[nsig], h[npars][npars],
		tmpmat[nsig], dZdif[nsig], Ndsig2[nsig], s1, s2, s3, s4, s5;		

	for(j = 0; j < npars; j++){
		for(i = 0; i < npars; i++){
			if(i <= j){						
				for(k = 0; k < npars; k++){
					derv1[k] = 0.0;
					derv2[k] = 0.0;
				}
				derv1[i] = 1.0;
				derv2[j] = 1.0;
				for(k = 0; k < nfact; k++){
					du1[k] = derv1[k];
					du2[k] = derv2[k];
				}
				for(k = nfact; k < npars; k++){
					dsig1[k-nfact] = derv1[k];
					dsig2[k-nfact] = derv2[k];
				}
				symMat(dsig1, &nfact);
				symMat(dsig2, &nfact);
				matrixMult(tmpmat, invSig, dsig2, &nfact); 
				matrixMult(dinvSig2, tmpmat, invSig, &nfact);				
				for(k = 0; k < nsig; k++)				
					dinvSig2[k] = -1.0 * dinvSig2[k];									
				outer(dZ, cMeans, du2, &nfact);				
				for(k = 0; k < nsig; k++)
					Ndsig2[k] = N * dsig2[k];
				matrixSub(dZdif, dZ, Ndsig2, &nfact);				
				matrixMult4(tmpmat, dsig1, dinvSig2, Zdif, invSig, &nfact);				
				s1 = 0.5 * tr(tmpmat, &nfact);
				matrixMult4(tmpmat, dsig1, invSig, Zdif, dinvSig2, &nfact);
				s2 = 0.5 * tr(tmpmat, &nfact);
				matrixMult4(tmpmat, dsig1, invSig, dZdif, invSig, &nfact);
				s3 = 0.5 * tr(tmpmat, &nfact);
				s4 = inner(du1, dinvSig2, cMeans, &nfact);
				s5 = N * inner(du1, invSig, du2, &nfact);				
				h[i][j] = s1 + s2 + s3 + s4 - s5;
				h[j][i] = h[i][j]; 
			}
		}
	}

	SEXP Rreturn;
	double *Preturn;
	PROTECT(Rreturn = allocMatrix(REALSXP,npars,npars));
	Preturn = NUMERIC_POINTER(Rreturn);	
	k=0;
	for(j = 0; j < npars; j++){
		for(i = 0; i < npars; i++){
			Preturn[k] = h[i][j];
			k++;
		}
	}
	
	UNPROTECT(7);	
	return(Rreturn);
}

