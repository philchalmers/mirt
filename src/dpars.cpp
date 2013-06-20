#include"Misc.h"

static double sumPrior(const NumericVector &a, const NumericVector &b)
{
    double ret = 0.0;
    for(int i = 0; i < a.length(); i++)
        ret += a(i) * b(i);
    return(ret);
}

static NumericVector makeOffterm(const NumericMatrix &dat, const NumericVector &p, const NumericVector &aTheta, 
        const double &D, const int &cat)
{
    NumericVector ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); CAT++){
        if(CAT == cat) continue;
        for(int n = 0; n < ret.length(); n++)
            ret(n) += dat(n, CAT) * p(n) * aTheta(n) * D;
    }
    return(ret);
}

static NumericVector makeOffterm2(const NumericMatrix &dat, const NumericVector &p1, const NumericVector &p2, 
        const NumericVector &aTheta, const double &D, const int &cat)
{
    NumericVector ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); CAT++){
        if(CAT == cat) continue;
        for(int n = 0; n < ret.length(); n++)
            ret(n) += dat(n, CAT) * p1(n) * p2(n) * aTheta(n) * D;
    }
    return(ret);
}

RcppExport SEXP dparsNominal(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta, SEXP RD,
        SEXP RPrior, SEXP RP, SEXP Rnum, SEXP Rdat, SEXP Rnfact, SEXP Rncat,
        SEXP Rakind, SEXP Rdind, SEXP Rak2, SEXP RP2, SEXP RP3,
        SEXP RaTheta, SEXP RaTheta2, SEXP Rdat_num, SEXP Rnumsum, SEXP RnumakD, 
        SEXP Rnumak2D2, SEXP RnumakDTheta_numsum, SEXP RestHess) 
{		
    BEGIN_RCPP
    IntegerVector Pnfact(Rnfact), Pncat(Rncat), Pakind(Rakind), Pdind(Rdind), 
                  estHess(RestHess);
    const int nfact = Pnfact(0);
    const int ncat = Pncat(0); 
    const int akind = Pakind(0); 
    const int dind = Pdind(0);
    int i, j, k, n;
    NumericVector a(Ra), ak(Rak), d(Rd), PD(RD), Prior(RPrior), ak2(Rak2),
                  numsum(Rnumsum), numakD(RnumakD), numak2D2(Rnumak2D2), 
                  aTheta(RaTheta), aTheta2(RaTheta2), dL(nfact + ncat*2);
    NumericVector unitNvec(aTheta.length()); 
    unitNvec.fill(1.0);
    const double D = PD(0); 
    const double D2 = PD(0)*PD(0);
    NumericMatrix Theta(RTheta), P(RP), num(Rnum), P2(RP2), P3(RP3),
                  dat_num(Rdat_num), numakDTheta_numsum(RnumakDTheta_numsum), 
                  d2L(nfact + ncat*2, nfact + ncat*2), dat(Rdat);
    const int N = dat.nrow();
    NumericVector tmpvec(N), tmpvec2(N), offterm(N), offterm2(N);
    
    //grad
    for(j = 0; j < nfact; j++){
        tmpvec.fill(0.0);
        for(i = 0; i < ncat; i++){
            for(n = 0; n < N; n++){
                tmpvec(n) += dat_num(n,i)*(D*ak(i)*Theta(n,j)*P(n,i) - 
                        P(n,i)*numakDTheta_numsum(n,j))*numsum(n);
            }
        }
        dL(j) = sumPrior(tmpvec, Prior);               
    }
    for(i = 0; i < ncat; i++){ 
        offterm = makeOffterm(dat, P(_,i), aTheta, D, i);
        offterm2 = makeOffterm(dat, P(_,i), unitNvec, D, i);
        for(n = 0; n < N; n++){ 
            tmpvec(n) = dat_num(n,i)*(D*aTheta(n)*P(n,i) - P2(n,i)*D*aTheta(n))*numsum(n) - offterm(n);
            tmpvec2(n) = dat_num(n,i)*(D*P(n,i) - P2(n,i)*D)*numsum(n) - offterm2(n);
        }
        dL(akind + i) = sumPrior(tmpvec, Prior);
        dL(dind + i) = sumPrior(tmpvec2, Prior);
    }

    //hess
    //a's
    if(estHess(0)){
        for(j = 0; j < nfact; j++){
            for(k = 0; k < nfact; k++){
                if(j <= k){
                    tmpvec.fill(0.0);
                    for(i = 0; i < ncat; i++){
                        for(n = 0; n < N; n++){
                            tmpvec(n) += dat_num(n,i)*(D2*ak2(i)*Theta(n,j)*Theta(n,k)*P(n,i) - 
                                    D*ak(i)*Theta(n,j)*P(n,i)*numakDTheta_numsum(n,k) -     
                                    D*ak(i)*Theta(n,k)*P(n,i)*numakDTheta_numsum(n,j) + 
                                    2*P(n,i)*numakD(n)*Theta(n,j)*numakD(n)*Theta(n,k)/ (numsum(n)*numsum(n)) - 
                                    P(n,i)*numak2D2(n)*Theta(n,j)*Theta(n,k)/numsum(n)) * numsum(n) - 
                                dat_num(n,i)*(D*ak(i)*Theta(n,j)*P(n,i) - P(n,i)*numakDTheta_numsum(n,j)) *
                                numsum(n)*D*ak(i)*Theta(n,k) + 
                                dat_num(n,i)*(D*ak(i)*Theta(n,j)*P(n,i) - P(n,i)*numakDTheta_numsum(n,j)) *
                                numakD(n)*Theta(n,k);
                        }
                    }
                    d2L(j,k) = sumPrior(tmpvec, Prior);
                    d2L(k, j) = d2L(j,k);
                }
            }
        }
        //a's with ak and d
        for(j = 0; j < nfact; j++){
            for(k = 0; k < ncat; k++){
                tmpvec.fill(0.0);
                tmpvec2.fill(0.0);
                for(i = 0; i < ncat; i++){
                    for(n = 0; n < N; n++){
                        if(i == k){
                            tmpvec(n) += dat_num(n,i)*(D2*ak(i)*Theta(n,j)*aTheta(n)*P(n,i) - 
                                        D*aTheta(n)*P(n,i)*numakDTheta_numsum(n,j) + 
                                        D*Theta(n,j)*P(n,i) - 2*D2*ak(i)*Theta(n,j)*aTheta(n)*P2(n,i) + 
                                        2*D*aTheta(n)*P2(n,i)*numakDTheta_numsum(n,j) - 
                                        D*Theta(n,j)*P2(n,i))*numsum(n) - 
                                dat_num(n,i)*(D*aTheta(n)*P(n,i) - D*aTheta(n)*P2(n,i))*numsum(n)*D*ak(i)*Theta(n,j) + 
                                dat_num(n,i)*(D*aTheta(n)*P(n,i) - D*aTheta(n)*P2(n,i))*(numakD(n)*Theta(n,j));
                        tmpvec2(n) += dat_num(n,i)*(D2*ak(i)*Theta(n,j)*P(n,i) - 
                                                            2*D2*ak(i)*Theta(n,j)*P2(n,i) -
                                                            D*P(n,i)*numakDTheta_numsum(n,j) +                                                          
                                                            2*D*P2(n,i)*numakDTheta_numsum(n,j))*numsum(n) - 
                            dat_num(n,i)*(D*P(n,i) - D*P2(n,i))*numsum(n)*D*ak(i)*Theta(n,j) + 
                            dat_num(n,i)*(D*P(n,i) - D*P2(n,i))*(numakD(n)*Theta(n,j));         
                        } else {
                            tmpvec(n) += -dat(n,i)*D2*ak(k)*aTheta(n)*Theta(n,j)*P(n,k) + 
                                dat(n,i)*P(n,k)*D*aTheta(n)*numakDTheta_numsum(n,j) - 
                                dat(n,i)*P(n,k)*D*Theta(n,j);
                            tmpvec2(n) += -dat(n,i)*D2*ak(k)*Theta(n,j)*P(n,k) + 
                                dat(n,i)*P(n,k)*D*numakDTheta_numsum(n,j);
                        }
                    }
                    d2L(j, akind + k) = sumPrior(tmpvec, Prior);
                    d2L(akind + k, j) = d2L(j, akind + k);
                    d2L(j, dind + k) = sumPrior(tmpvec2, Prior);
                    d2L(dind + k, j) = d2L(j, dind + k);
                }
            }
        }
        //ak's and d's
        for(j = 0; j < ncat; j++){
            tmpvec = makeOffterm(dat, P2(_,j), aTheta2, D2, j);
            tmpvec2 = makeOffterm(dat, P(_,j), aTheta2, D2, j);
            for(n = 0; n < N; n++)
                offterm(n) = tmpvec(n) - tmpvec2(n);
            tmpvec = makeOffterm(dat, P2(_,j), unitNvec, D2, j);
            tmpvec2 = makeOffterm(dat, P(_,j), unitNvec, D2, j);
            for(n = 0; n < N; n++)
                offterm2(n) = tmpvec(n) - tmpvec2(n);
            for(n = 0; n < N; n++){
                tmpvec(n) = dat_num(n,j)*(D2*aTheta2(n)*P(n,j) - 3*D2*aTheta2(n)*P2(n,j) + 
                                                2*D2*aTheta2(n)*P3(n,j))*numsum(n) - dat(n,j)/num(n,j)*(D*aTheta(n)*P(n,j) - 
                                                D*aTheta(n)*P2(n,j))*numsum(n)*D*aTheta(n) + dat(n,j)*(D*aTheta(n)*P(n,j) - 
                                                D*aTheta(n)*P2(n,j))*D*aTheta(n) + offterm(n);
                tmpvec2(n) = dat_num(n,j)*(D2*P(n,j) - 3*D2*P2(n,j) + 
                                            2*D2*P3(n,j))*numsum(n) - dat(n,j)/num(n,j)*(D*P(n,j) - 
                                            D*P2(n,j))*numsum(n)*D + dat(n,j)*(D*P(n,j) - 
                                            D*P2(n,j))*D + offterm2(n);
            }
            d2L(akind + j, akind + j) = sumPrior(tmpvec, Prior); 
            d2L(dind + j, dind + j) = sumPrior(tmpvec2, Prior);
            for(i = 0; i < ncat; i++){
                if(j < i){   
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta2, D2, i);
                    offterm2 = makeOffterm2(dat, P(_,j), P(_,i), unitNvec, D2, i);
                    for(n = 0; n < N; n++){
                        tmpvec(n) = dat_num(n,i) * (-D2*aTheta2(n)*P(n,i)*P(n,j) + 2*P2(n,i) *D2*aTheta2(n)*P(n,j))*numsum(n) + 
                                     dat_num(n,i) * (D*aTheta(n)*P(n,i) - P2(n,i) * D * aTheta(n))*D*aTheta(n)*num(n,j)+offterm(n);
                        tmpvec2(n) = dat_num(n,i) * (-D2*P(n,i)*P(n,j) + 2*P2(n,i) *D2*P(n,j)) * numsum(n) + 
                            dat_num(n,i) * (D*P(n,i) - P2(n,i)*D)*D * num(n,j) + offterm2(n);
                    }
                    d2L(akind + i, akind + j) = sumPrior(tmpvec, Prior);
                    d2L(akind + j, akind + i) = d2L(akind + i, akind + j);
                    d2L(dind + i, dind + j) = sumPrior(tmpvec2, Prior); 
                    d2L(dind + j, dind + i) = d2L(dind + i, dind + j);
                }
                if(abs(j-i) == 0){
                    tmpvec = makeOffterm(dat, P2(_,i), aTheta, D2, i);
                    tmpvec2 = makeOffterm(dat, P(_,i), aTheta, D2, i);
                    for(n = 0; n < N; n++){
                        offterm(n) = tmpvec(n) - tmpvec2(n);
                        tmpvec(n) = dat_num(n,i)*(D2*aTheta(n)*P(n,i) - 3*D2*aTheta(n)*P2(n,i) + 
                                2*D2*aTheta(n)*P3(n,i))*numsum(n) - dat_num(n,i)*(D*aTheta(n)*P(n,i) - 
                                D*aTheta(n)*P2(n,i))*numsum(n)*D + dat(n,i)*(D*P(n,i) - 
                                D*P2(n,i))*D*aTheta(n) + offterm(n);
                    }
                    d2L(dind + j, akind + i) = sumPrior(tmpvec, Prior);
                    d2L(akind + i, dind + j) = d2L(dind + j, akind + i);
                } else {
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta, D2, i);
                    for(n = 0; n < N; n++){
                        tmpvec(n) = dat_num(n,i) * (-D2*aTheta(n)*P(n,i)*P(n,j) + 2*P2(n,i) *D2*aTheta(n)*P(n,j)) * numsum(n) + 
                            dat_num(n,i) * (D*P(n,i) - P2(n,i) * D) * D * aTheta(n) * num(n,j) + offterm(n);                
                    }
                    d2L(akind + i, dind + j) = sumPrior(tmpvec, Prior);
                    d2L(dind + j, akind + i) = d2L(akind + i, dind + j);
                }
            }
        }            
    }

    List ret;
    ret["grad"] = dL;
    ret["hess"] = d2L;
	return(ret);
	END_RCPP
}


RcppExport SEXP dparsPoly(SEXP Rprob, SEXP RThetas, SEXP RPrior, SEXP Rdat, SEXP Rnzeta, SEXP RestHess) 
{		
    BEGIN_RCPP
    /* 
        Rprob = numeric matrix of probabilities
        RThetas = numeric matrix of abilities
        Rdat = integer matrix of dichotomized item responses
        nzeta = integer number of response categories
     */

	int i, j, k; 
	NumericMatrix prob(Rprob);
	NumericMatrix Thetas(RThetas);
    NumericVector Prior(RPrior);
    NumericMatrix dat2(Rdat);
    IntegerVector Pnzeta(Rnzeta);
    IntegerVector estHess(RestHess);
    const int nzeta = Pnzeta[0];
    const int nfact = Thetas.ncol();
    const int N = Thetas.nrow();
    NumericMatrix dat(dat2.nrow(), dat2.ncol());
    NumericMatrix d2L(nfact + nzeta, nfact + nzeta);
    NumericVector dL(nfact + nzeta);

	NumericVector Pk(N), Pk_1(N), Pk_p1(N), PQ_1(N), PQ(N), PQ_p1(N), 
			Pk_1Pk(N), Pk_Pkp1(N), dif1(N), dif1sq(N), dif2(N), 
			dif2sq(N), tmp1(N), tmp2(N), tmp3(N), csums(nfact);			
	NumericMatrix P(N,nzeta+2), PQfull(N,nzeta+2), mattmp(N,nfact), d2Louter;	   
	double tmp;
	IntegerVector factind(nfact);
	for(j = 0; j < (nzeta + 2); j++){
		for(i = 0; i < N; i++){
			P(i,j) = prob(i,j);
			PQfull(i,j) = prob(i,j) * (1.0 - prob(i,j));
		}
	}
	for(j = 0; j < dat2.ncol(); j++){
		for(i = 0; i < N; i++){
		    dat(i,j) = dat2(i,j) * Prior(i);
		}
	}
	for(j = 0; j < nfact; j++)
		factind(j) = nzeta + j;
	for(j = 0; j < (nzeta + 1); j++){
		if(j < nzeta){
			for(i = 0; i < N; i++){
				Pk_1(i) = P(i,j);
				Pk(i) = P(i,j + 1);
				Pk_p1(i) = P(i,j + 2);
				PQ_1(i) = PQfull(i,j);
				PQ(i) = PQfull(i,j + 1);
				PQ_p1(i) = PQfull(i,j + 2);
				Pk_1Pk(i) = Pk_1(i) - Pk(i);
				Pk_Pkp1(i) = Pk(i) - Pk_p1(i);
                if(Pk_1Pk(i) < 1e-10) Pk_1Pk(i) = 1e-10;
                if(Pk_Pkp1(i) < 1e-10) Pk_Pkp1(i) = 1e-10;
				dif1(i) = dat(i,j) / Pk_1Pk(i);
				dif1sq(i) = dat(i,j) / (Pk_1Pk(i) * Pk_1Pk(i));
				dif2(i) = dat(i,j+1) / Pk_Pkp1(i);
				dif2sq(i) = dat(i,j+1) / (Pk_Pkp1(i) * Pk_Pkp1(i));
			}			
			tmp = 0.0;
			for(i = 0; i < N; i++)
				tmp += (-1.0) * PQ(i) * (dif1(i) - dif2(i));			
			dL(j) = tmp;			
			if(estHess(0)){
			    tmp = 0.0;
			    for(i = 0; i < N; i++)
			    	tmp += (-1.0) * PQ(i) * PQ(i) * (dif1sq(i) + dif2sq(i)) -				
			    		(dif1(i) - dif2(i)) * (Pk(i) * (1.0 - Pk(i)) * (1.0 - 2.0*Pk(i)));			
			    d2L(j,j) = tmp;
			    if(j < (nzeta - 1)){
			    	tmp = 0.0;
			    	for(i = 0; i < N; i++)
			    		tmp += dif2sq(i) * PQ_p1(i) * PQ(i);
			    	d2L(j,j+1) = tmp;
			    	d2L(j+1,j) = tmp;
			    }
			    for(i = 0; i < N; i++){
			    	tmp1(i) = (-1.0) * dif2sq(i) * PQ(i) * (PQ(i) - PQ_p1(i));
			    	tmp2(i) = dif1sq(i) * PQ(i) * (PQ_1(i) - PQ(i));
			    	tmp3(i) = (dif1(i) - dif2(i)) * (Pk(i) * (1.0 - Pk(i)) * (1.0 - 2.0*Pk(i)));
			    }
			    for(k = 0; k < nfact; k++){
			    	csums(k) = 0.0;
			    	for(i = 0; i < N; i++){
			    		mattmp(i,k) = tmp1(i) * Thetas(i,k) + tmp2(i) * Thetas(i,k) - 
			    			tmp3(i) * Thetas(i,k);
			    		csums(k) += mattmp(i,k);
			    	}
			    }
			    for(i = 0; i < nfact; i++){
			    	d2L(j,factind(i)) = csums(i);
			    	d2L(factind(i),j) = csums(i);
			    }			
            }
		} else {					
			for(i = 0; i < N; i++){
				Pk_1(i) = P(i,j);
				Pk(i) = P(i,j + 1);			
				PQ_1(i) = PQfull(i,j);
				PQ(i) = PQfull(i,j + 1);			
				Pk_1Pk(i) = Pk_1(i) - Pk(i);
                if(Pk_1Pk(i) < 1e-10) Pk_1Pk(i) = 1e-10; 
				dif1(i) = dat(i,j) / Pk_1Pk(i);
				dif1sq(i) = dat(i,j) / (Pk_1Pk(i) * Pk_1Pk(i));			
			}	
		}
		for(k = 0; k < nfact; k++){
			csums(k) = 0.0;
			for(i = 0; i < N; i++){
				mattmp(i,k) = dif1(i) * (PQ_1(i) - PQ(i)) * Thetas(i,k);
				csums(k) += mattmp(i,k);
			}
		}
		for(i = 0; i < nfact; i++)
    		dL(factind(i)) += csums(i);			
		
		if(estHess(0)){
		    d2Louter = polyOuter(Thetas, Pk, Pk_1, PQ_1, PQ, dif1sq, dif1);		
		    for(k = 0; k < nfact; k++)
			    for(i = 0; i < nfact; i++)
				    d2L(factind(i),factind(k)) += d2Louter(i,k);				
        }
	}

    List ret;
    ret["grad"] = dL;
    ret["hess"] = d2L;
	return(ret);
	END_RCPP
}

RcppExport SEXP dparsDich(SEXP Ra, SEXP Rd, SEXP Rg, SEXP Ru, SEXP RD, SEXP RTheta,
    SEXP RPrior, SEXP Rr1, SEXP Rr2, SEXP RestHess, SEXP asMatrix, SEXP Rot) 
{		
    BEGIN_RCPP
    
    int i, j;     
    NumericVector a(Ra);
    NumericVector d(Rd);
    NumericVector Pg(Rg);
    NumericVector Pu(Ru);
    NumericVector PD(RD);
	NumericMatrix Theta(RTheta);
    NumericVector Prior(RPrior);
    NumericVector r1(Rr1);
    NumericVector r2(Rr2);    
    NumericVector ot(Rot);
    IntegerVector estHess(RestHess);        
    const int nfact = Theta.ncol();    
    NumericVector P, Pstar, Q, Qstar;
    const double g = Pg(0);
    const double u = Pu(0);
    const double D = PD(0);
    const double g0 = 0.0;
    const double u1 = 1.0;    
    NumericMatrix hess (nfact + 3, nfact + 3);
    NumericVector grad (nfact + 3);
    NumericVector r1_P, r1_P2, r2_Q2, r2_Q;    
    
    P = itemTrace(a, &d(0), Theta, &g, &u, &D, ot);
    Pstar = itemTrace(a, &d(0), Theta, &g0, &u1, &D, ot);
    Q = 1.0 - P;
    Qstar = 1.0 - Pstar;        
    r1_P = r1/P;
    r1_P2 = r1/(P*P);
    r2_Q = r2/Q; 
    r2_Q2 = r2/(Q*Q);
    grad(nfact) = sum((u-g)*D*Pstar*Qstar*(r1_P - r2_Q)*Prior);    
    grad(nfact + 1) = sum(Qstar*(r1_P - r2_Q)*Prior);
    grad(nfact + 2) = sum(Pstar*(r1_P - r2_Q)*Prior);
    for(i = 0; i < nfact; i++)
        grad(i) = sum(Theta(_, i) * D * Pstar * Qstar * (u-g) * (r1_P - r2_Q) * Prior);
        
    if(estHess){
        int gloc = nfact+1; 
        int uloc = nfact+2;
        double ugD2 = (u-g) * D*D;
        double ugD = (u-g) * D; 
        NumericVector Pstar2 = Pstar*Pstar; 
        NumericVector Pstar3 = Pstar*Pstar*Pstar;
        hess(nfact,nfact) = sum((r1_P * (ugD2 * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                          r1_P2 * (ugD * (Pstar - Pstar2))*(ugD * (Pstar - Pstar2)) +
                                          r2_Q * (ugD2 * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                          r2_Q2 * (ugD * (-Pstar + Pstar2))*(ugD * (-Pstar + Pstar2))) * Prior);
        hess(gloc,gloc) = -1.0 * sum(Qstar*Qstar *(r1_P2 + r2_Q2)*Prior);
        hess(uloc,uloc) = -1.0 * sum(Pstar2 *(r1_P2 + r2_Q2)*Prior);
        hess(nfact, gloc) = sum((r1_P * (D * (-Pstar + Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Qstar +
                     r2_Q * (D * (Pstar - Pstar2)) -
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * -Qstar) * Prior);
        hess(gloc, nfact) = hess(nfact, gloc);
        hess(nfact, uloc) = sum((r1_P * (D * (Pstar - Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Pstar +
                     r2_Q * (D * (-Pstar + Pstar2)) +
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * Pstar) * Prior);
        hess(uloc, nfact) = hess(nfact, uloc);
        hess(gloc, uloc) = sum((-1.0 * r1_P2 * Pstar * Qstar + r2_Q2 * Pstar * (-1.0 + Pstar )) * Prior);
        hess(uloc, gloc) = hess(gloc, uloc);
        for(i = 0; i < nfact; i++)
            for(j = 0; j < nfact; j++)
                if(i <= j)
                    hess(i, j) = sum((r1_P * (ugD2 * Theta(_,i) * Theta(_,j) * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                           r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) *
                                              (ugD * Theta(_,j) * (Pstar - Pstar2)) +
                                           r2_Q * (ugD2 * Theta(_,i) * Theta(_,j) *
                                               (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                           r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2)) *
                                              (ugD * Theta(_,j) * (-Pstar + Pstar2))) * Prior);                    
        for(i = 0; i < nfact; i++){
            hess(i, nfact) = 
                sum((r1_P * (ugD2 * Theta(_,i) * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) *
                            (ugD * (Pstar - Pstar2)) +
                         r2_Q * (ugD2 * Theta(_,i) * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2)) *
                         (ugD * (-Pstar + Pstar2))) * Prior);
            hess(nfact, i) = hess(i, nfact);
            hess(i, gloc) = 
                sum((r1_P * (D * Theta(_,i) * (-Pstar + Pstar2)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) * Qstar +
                         r2_Q * (D * Theta(_,i) * (Pstar - Pstar2)) -
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2) ) * (Pstar - 1.0)) * Prior);
            hess(gloc, i) = hess(i, gloc);
            hess(i, uloc) =  
                sum((r1_P * (D * Theta(_,i) * (Pstar - Pstar2)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) * Pstar +
                         r2_Q * (D * Theta(_,i) * (-Pstar + Pstar2)) +
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2) ) * Pstar) * Prior);
            hess(uloc, i) = hess(i, uloc);
        }
    }    
    List ret;    
    ret["grad"] = grad;
    ret["hess"] = hess;
	return(ret);
	END_RCPP
}
