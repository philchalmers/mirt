#include"Misc.h"

static NumericVector makeOffterm(const NumericMatrix &dat, const NumericVector &p, const NumericVector &aTheta, 
        const int &cat)
{
    NumericVector ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); CAT++){
        if(CAT == cat) continue;
        for(int n = 0; n < ret.length(); ++n)
            ret(n) += dat(n, CAT) * p(n) * aTheta(n);
    }
    return(ret);
}

static NumericVector makeOffterm2(const NumericMatrix &dat, const NumericVector &p1, const NumericVector &p2, 
        const NumericVector &aTheta, const int &cat)
{
    NumericVector ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); CAT++){
        if(CAT == cat) continue;
        for(int n = 0; n < ret.length(); ++n)
            ret(n) += dat(n, CAT) * p1(n) * p2(n) * aTheta(n);
    }
    return(ret);
}

static double difexp(const double *x)
{
    double ret = *x * (1.0 - *x);
    return(ret);
}

RcppExport SEXP dparsNominal(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta,
        SEXP RP, SEXP Rnum, SEXP Rdat, SEXP Rnfact, SEXP Rncat,
        SEXP Rakind, SEXP Rdind, SEXP Rak2, SEXP RP2, SEXP RP3,
        SEXP RaTheta, SEXP RaTheta2, SEXP Rdat_num, SEXP Rnumsum, SEXP RnumakD, 
        SEXP Rnumak2D2, SEXP RnumakDTheta_numsum, SEXP RestHess) 
{		
    BEGIN_RCPP
    const int nfact = as<int>(Rnfact);
    const int ncat = as<int>(Rncat); 
    const int akind = as<int>(Rakind);
    const int dind = as<int>(Rdind);
    const int estHess = as<int>(RestHess);
    const NumericVector a(Ra), ak(Rak), d(Rd), ak2(Rak2),
                  numsum(Rnumsum), numakD(RnumakD), numak2D2(Rnumak2D2), 
                  aTheta(RaTheta), aTheta2(RaTheta2);
    NumericVector dL(nfact + ncat*2);
    const NumericVector unitNvec(aTheta.length(), 1.0);
    const NumericMatrix Theta(RTheta), P(RP), num(Rnum), P2(RP2), P3(RP3),
                  dat_num(Rdat_num), numakDTheta_numsum(RnumakDTheta_numsum), dat(Rdat);
    NumericMatrix d2L(nfact + ncat*2, nfact + ncat*2);
    const int N = dat.nrow();
    NumericVector tmpvec(N), tmpvec2(N), offterm(N), offterm2(N);
    
    //grad
    for(int j = 0; j < nfact; ++j){
        tmpvec.fill(0.0);
        for(int i = 0; i < ncat; ++i){
            for(int n = 0; n < N; ++n){
                tmpvec(n) += dat_num(n,i)*(ak(i)*Theta(n,j)*P(n,i) - 
                        P(n,i)*numakDTheta_numsum(n,j))*numsum(n);
            }
        }
        dL(j) = sum(tmpvec);               
    }
    for(int i = 0; i < ncat; ++i){ 
        offterm = makeOffterm(dat, P(_,i), aTheta,  i);
        offterm2 = makeOffterm(dat, P(_,i), unitNvec,  i);
        for(int n = 0; n < N; ++n){ 
            tmpvec(n) = dat_num(n,i)*(aTheta(n)*P(n,i) - P2(n,i)*aTheta(n))*numsum(n) - offterm(n);
            tmpvec2(n) = dat_num(n,i)*(P(n,i) - P2(n,i))*numsum(n) - offterm2(n);
        }
        dL(akind + i) = sum(tmpvec);
        dL(dind + i) = sum(tmpvec2);
    }

    //hess
    //a's
    if(estHess){
        for(int j = 0; j < nfact; ++j){
            for(int k = 0; k < nfact; ++k){
                if(j <= k){
                    tmpvec.fill(0.0);
                    for(int i = 0; i < ncat; ++i){
                        for(int n = 0; n < N; ++n){
                            tmpvec(n) += dat_num(n,i)*(ak2(i)*Theta(n,j)*Theta(n,k)*P(n,i) - 
                                    ak(i)*Theta(n,j)*P(n,i)*numakDTheta_numsum(n,k) -     
                                    ak(i)*Theta(n,k)*P(n,i)*numakDTheta_numsum(n,j) + 
                                    2*P(n,i)*numakD(n)*Theta(n,j)*numakD(n)*Theta(n,k)/ (numsum(n)*numsum(n)) - 
                                    P(n,i)*numak2D2(n)*Theta(n,j)*Theta(n,k)/numsum(n)) * numsum(n) - 
                                dat_num(n,i)*(ak(i)*Theta(n,j)*P(n,i) - P(n,i)*numakDTheta_numsum(n,j)) *
                                numsum(n)*ak(i)*Theta(n,k) + 
                                dat_num(n,i)*(ak(i)*Theta(n,j)*P(n,i) - P(n,i)*numakDTheta_numsum(n,j)) *
                                numakD(n)*Theta(n,k);
                        }
                    }
                    d2L(j,k) = sum(tmpvec);
                    d2L(k, j) = d2L(j,k);
                }
            }
        }
        //a's with ak and d
        for(int j = 0; j < nfact; ++j){
            for(int k = 0; k < ncat; ++k){
                tmpvec.fill(0.0);
                tmpvec2.fill(0.0);
                for(int i = 0; i < ncat; ++i){
                    for(int n = 0; n < N; ++n){
                        if(i == k){
                            tmpvec(n) += dat_num(n,i)*(ak(i)*Theta(n,j)*aTheta(n)*P(n,i) - 
                                        aTheta(n)*P(n,i)*numakDTheta_numsum(n,j) + 
                                        Theta(n,j)*P(n,i) - 2*ak(i)*Theta(n,j)*aTheta(n)*P2(n,i) + 
                                        2*aTheta(n)*P2(n,i)*numakDTheta_numsum(n,j) - 
                                        Theta(n,j)*P2(n,i))*numsum(n) - 
                                dat_num(n,i)*(aTheta(n)*P(n,i) - aTheta(n)*P2(n,i))*numsum(n)*ak(i)*Theta(n,j) + 
                                dat_num(n,i)*(aTheta(n)*P(n,i) - aTheta(n)*P2(n,i))*(numakD(n)*Theta(n,j));
                        tmpvec2(n) += dat_num(n,i)*(ak(i)*Theta(n,j)*P(n,i) - 
                                                            2*ak(i)*Theta(n,j)*P2(n,i) -
                                                            P(n,i)*numakDTheta_numsum(n,j) +                                                          
                                                            2*P2(n,i)*numakDTheta_numsum(n,j))*numsum(n) - 
                            dat_num(n,i)*(P(n,i) - P2(n,i))*numsum(n)*ak(i)*Theta(n,j) + 
                            dat_num(n,i)*(P(n,i) - P2(n,i))*(numakD(n)*Theta(n,j));         
                        } else {
                            tmpvec(n) += -dat(n,i)*ak(k)*aTheta(n)*Theta(n,j)*P(n,k) + 
                                dat(n,i)*P(n,k)*aTheta(n)*numakDTheta_numsum(n,j) - 
                                dat(n,i)*P(n,k)*Theta(n,j);
                            tmpvec2(n) += -dat(n,i)*ak(k)*Theta(n,j)*P(n,k) + 
                                dat(n,i)*P(n,k)*numakDTheta_numsum(n,j);
                        }
                    }
                    d2L(j, akind + k) = sum(tmpvec);
                    d2L(akind + k, j) = d2L(j, akind + k);
                    d2L(j, dind + k) = sum(tmpvec2);
                    d2L(dind + k, j) = d2L(j, dind + k);
                }
            }
        }
        //ak's and d's
        for(int j = 0; j < ncat; ++j){
            tmpvec = makeOffterm(dat, P2(_,j), aTheta2, j);
            tmpvec2 = makeOffterm(dat, P(_,j), aTheta2, j);
            for(int n = 0; n < N; ++n)
                offterm(n) = tmpvec(n) - tmpvec2(n);
            tmpvec = makeOffterm(dat, P2(_,j), unitNvec, j);
            tmpvec2 = makeOffterm(dat, P(_,j), unitNvec, j);
            for(int n = 0; n < N; ++n)
                offterm2(n) = tmpvec(n) - tmpvec2(n);
            for(int n = 0; n < N; ++n){
                tmpvec(n) = dat_num(n,j)*(aTheta2(n)*P(n,j) - 3*aTheta2(n)*P2(n,j) + 
                                                2*aTheta2(n)*P3(n,j))*numsum(n) - dat(n,j)/num(n,j)*(aTheta(n)*P(n,j) - 
                                                aTheta(n)*P2(n,j))*numsum(n)*aTheta(n) + dat(n,j)*(aTheta(n)*P(n,j) - 
                                                aTheta(n)*P2(n,j))*aTheta(n) + offterm(n);
                tmpvec2(n) = dat_num(n,j)*(P(n,j) - 3*P2(n,j) + 
                                            2*P3(n,j))*numsum(n) - dat(n,j)/num(n,j)*(P(n,j) - 
                                            P2(n,j))*numsum(n) + dat(n,j)*(P(n,j) - 
                                            P2(n,j)) + offterm2(n);
            }
            d2L(akind + j, akind + j) = sum(tmpvec); 
            d2L(dind + j, dind + j) = sum(tmpvec2);
            for(int i = 0; i < ncat; ++i){
                if(j < i){   
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta2, i);
                    offterm2 = makeOffterm2(dat, P(_,j), P(_,i), unitNvec, i);
                    for(int n = 0; n < N; ++n){
                        tmpvec(n) = dat_num(n,i) * (-aTheta2(n)*P(n,i)*P(n,j) + 2*P2(n,i) *aTheta2(n)*P(n,j))*numsum(n) + 
                                     dat_num(n,i) * (aTheta(n)*P(n,i) - P2(n,i) * aTheta(n))*aTheta(n)*num(n,j)+offterm(n);
                        tmpvec2(n) = dat_num(n,i) * (-P(n,i)*P(n,j) + 2*P2(n,i) *P(n,j)) * numsum(n) + 
                            dat_num(n,i) * (P(n,i) - P2(n,i)) * num(n,j) + offterm2(n);
                    }
                    d2L(akind + i, akind + j) = sum(tmpvec);
                    d2L(akind + j, akind + i) = d2L(akind + i, akind + j);
                    d2L(dind + i, dind + j) = sum(tmpvec2); 
                    d2L(dind + j, dind + i) = d2L(dind + i, dind + j);
                }
                if(abs(j-i) == 0){
                    tmpvec = makeOffterm(dat, P2(_,i), aTheta, i);
                    tmpvec2 = makeOffterm(dat, P(_,i), aTheta, i);
                    for(int n = 0; n < N; ++n){
                        offterm(n) = tmpvec(n) - tmpvec2(n);
                        tmpvec(n) = dat_num(n,i)*(aTheta(n)*P(n,i) - 3*aTheta(n)*P2(n,i) + 
                                2*aTheta(n)*P3(n,i))*numsum(n) - dat_num(n,i)*(aTheta(n)*P(n,i) - 
                                aTheta(n)*P2(n,i))*numsum(n) + dat(n,i)*(P(n,i) - 
                                P2(n,i))*aTheta(n) + offterm(n);
                    }
                    d2L(dind + j, akind + i) = sum(tmpvec);
                    d2L(akind + i, dind + j) = d2L(dind + j, akind + i);
                } else {
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta, i);
                    for(int n = 0; n < N; ++n){
                        tmpvec(n) = dat_num(n,i) * (-aTheta(n)*P(n,i)*P(n,j) + 2*P2(n,i) *aTheta(n)*P(n,j)) * numsum(n) + 
                            dat_num(n,i) * (P(n,i) - P2(n,i)) * aTheta(n) * num(n,j) + offterm(n);                
                    }
                    d2L(akind + i, dind + j) = sum(tmpvec);
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


RcppExport SEXP dparsPoly(SEXP Rprob, SEXP RThetas, SEXP Rdat, SEXP Rnzeta, SEXP RestHess) 
{		
    BEGIN_RCPP    

	const NumericMatrix prob(Rprob);
	const NumericMatrix Thetas(RThetas);    
    const NumericMatrix dat(Rdat);
    const int nzeta = as<int>(Rnzeta);
    const int estHess = as<int>(RestHess);        
    const int nfact = Thetas.ncol();
    const int N = Thetas.nrow();
    NumericMatrix d2L(nfact + nzeta, nfact + nzeta);
    NumericVector dL(nfact + nzeta);

	vector<double> Pk(N), Pk_1(N), Pk_p1(N), PQ_1(N), PQ(N), PQ_p1(N), 
			Pk_1Pk(N), Pk_Pkp1(N), dif1(N), dif1sq(N), dif2(N), 
			dif2sq(N), tmp1(N), tmp2(N), tmp3(N), csums(nfact);
	NumericMatrix P(N,nzeta+2), PQfull(N,nzeta+2), mattmp(N,nfact), d2Louter;
	IntegerVector factind(nfact);
	for(int j = 0; j < (nzeta + 2); ++j){
		for(int i = 0; i < N; ++i){
			P(i,j) = prob(i,j);
			PQfull(i,j) = prob(i,j) * (1.0 - prob(i,j));
		}
	}	
	for(int j = 0; j < nfact; ++j)
		factind(j) = nzeta + j;
	for(int j = 0; j < (nzeta + 1); ++j){
		if(j < nzeta){
			for(int i = 0; i < N; ++i){
				Pk_1[i] = P(i,j);
				Pk[i] = P(i,j + 1);
				Pk_p1[i] = P(i,j + 2);
				PQ_1[i] = PQfull(i,j);
				PQ[i] = PQfull(i,j + 1);
				PQ_p1[i] = PQfull(i,j + 2);
				Pk_1Pk[i] = Pk_1[i] - Pk[i];
				Pk_Pkp1[i] = Pk[i] - Pk_p1[i];
                if(Pk_1Pk[i] < 1e-10) Pk_1Pk[i] = 1e-10;
                if(Pk_Pkp1[i] < 1e-10) Pk_Pkp1[i] = 1e-10;
				dif1[i] = dat(i,j) / Pk_1Pk[i];
				dif1sq[i] = dat(i,j) / (Pk_1Pk[i] * Pk_1Pk[i]);
				dif2[i] = dat(i,j+1) / Pk_Pkp1[i];
				dif2sq[i] = dat(i,j+1) / (Pk_Pkp1[i] * Pk_Pkp1[i]);
			}		
			for(int i = 0; i < N; ++i)
				dL(j) += (-1.0) * PQ[i] * (dif1[i] - dif2[i]);			
			if(estHess){
			    for(int i = 0; i < N; ++i)
			    	d2L(j,j) += (-1.0) * PQ[i] * PQ[i] * (dif1sq[i] + dif2sq[i]) -				
			    		(dif1[i] - dif2[i]) * (Pk[i] * (1.0 - Pk[i]) * (1.0 - 2.0*Pk[i]));			
			    if(j < (nzeta - 1)){
			    	for(int i = 0; i < N; ++i)
			    		d2L(j,j+1) += dif2sq[i] * PQ_p1[i] * PQ[i];
			    	d2L(j+1,j) = d2L(j,j+1);
			    }
			    for(int i = 0; i < N; ++i){
			    	tmp1[i] = (-1.0) * dif2sq[i] * PQ[i] * (PQ[i] - PQ_p1[i]);
			    	tmp2[i] = dif1sq[i] * PQ[i] * (PQ_1[i] - PQ[i]);
			    	tmp3[i] = (dif1[i] - dif2[i]) * (Pk[i] * (1.0 - Pk[i]) * (1.0 - 2.0*Pk[i]));
			    }
			    for(int k = 0; k < nfact; ++k){
			    	csums[k] = 0.0;
			    	for(int i = 0; i < N; ++i){
			    		mattmp(i,k) = tmp1[i] * Thetas(i,k) + tmp2[i] * Thetas(i,k) - 
			    			tmp3[i] * Thetas(i,k);
			    		csums[k] += mattmp(i,k);
			    	}
			    }
			    for(int i = 0; i < nfact; ++i){
			    	d2L(j,factind(i)) = csums[i];
			    	d2L(factind(i),j) = csums[i];
			    }			
            }
		} else {					
			for(int i = 0; i < N; ++i){
				Pk_1[i] = P(i,j);
				Pk[i] = P(i,j + 1);			
				PQ_1[i] = PQfull(i,j);
				PQ[i] = PQfull(i,j + 1);			
				Pk_1Pk[i] = Pk_1[i] - Pk[i];
                if(Pk_1Pk[i] < 1e-10) Pk_1Pk[i] = 1e-10; 
				dif1[i] = dat(i,j) / Pk_1Pk[i];
				dif1sq[i] = dat(i,j) / (Pk_1Pk[i] * Pk_1Pk[i]);			
			}	
		}
		for(int k = 0; k < nfact; ++k){
			csums[k] = 0.0;
			for(int i = 0; i < N; ++i){
				mattmp(i,k) = dif1[i] * (PQ_1[i] - PQ[i]) * Thetas(i,k);
				csums[k] += mattmp(i,k);
			}
		}
		for(int i = 0; i < nfact; ++i)
    		dL(factind(i)) += csums[i];			
		
		if(estHess){
		    d2Louter = polyOuter(Thetas, Pk, Pk_1, PQ_1, PQ, dif1sq, dif1);		
		    for(int k = 0; k < nfact; ++k)
			    for(int i = 0; i < nfact; ++i)
				    d2L(factind(i),factind(k)) += d2Louter(i,k);				
        }
	}
    
    //reorder 
    NumericVector grad(dL.length());
    NumericMatrix hess(d2L.ncol(), d2L.ncol());
    
    for(int i = 0; i < nfact; ++i)
        grad(i) = dL(i+nzeta);
    for(int i = 0; i < nzeta; ++i)
        grad(i+nfact) = dL(i);    
    if(estHess){
        for(int i = 0; i < nfact; ++i){
            for(int j = 0; j < nfact; ++j){
                hess(i,j) = d2L(i+nzeta, j+nzeta);
            }
        }
        for(int i = 0; i < nzeta; ++i){
            for(int j = 0; j < nzeta; ++j){
                hess(i+nfact,j+nfact) = d2L(i, j);
            }
        }
        for(int i = 0; i < nfact; ++i){
            for(int j = 0; j < nzeta; ++j){
                hess(j+nfact, i) = d2L(nzeta+i,j);
                hess(i, j+nfact) = d2L(nzeta+i,j);
            }
        }        
    }

    List ret;
    ret["grad"] = grad;
    ret["hess"] = hess;      
	return(ret);
	END_RCPP
}

RcppExport SEXP dparsDich(SEXP Rx, SEXP RTheta, SEXP RestHess, SEXP Rrs, SEXP Rot) 
{		
    BEGIN_RCPP
    
	const NumericVector par(Rx);
    const NumericMatrix Theta(RTheta);
	const NumericMatrix rs(Rrs);	
	const NumericVector ot(Rot);
    const int estHess = as<int>(RestHess);
			
	const int nfact = Theta.ncol();    
    const int N = Theta.nrow();
	vector<double> a(nfact);
	for(int i = 0; i < nfact; ++i) a[i] = par(i);
	const double d = par(nfact);
    const double expg = par(nfact+1);
    const double expu = par(nfact+2);    
    const double g = antilogit(&expg);
    const double u = antilogit(&expu);
    const double difexpg = difexp(&g);
    const double difexpu = difexp(&u);
    const double ugD = (u-g); 
    const double gm1 = (1.0 - g);
    const double um1 = (1.0 - u);
    const double u_1u = u * um1;
    const double g_1g = g * gm1;
    const int gloc = nfact+1; 
    const int uloc = nfact+2;
    NumericMatrix hess (nfact + 3, nfact + 3);
    vector<double> grad (nfact + 3);
    
    vector<double> P(N), Pstar(N);
    itemTrace(P, Pstar, a, &d, Theta, &g, &u, ot);
    
    for(int i = 0; i < N; ++i){
        double Q = 1.0 - P[i];
        double Qstar = 1.0 - Pstar[i];
        double r1_P = rs(i, 1) / P[i];
        double r1_P2 = rs(i, 1) / (P[i]*P[i]);
        double r2_Q = rs(i, 0) / Q; 
        double r2_Q2 = rs(i, 0) / (Q*Q);
        double r1_Pr2_Q = r1_P - r2_Q;
        grad[nfact] += (u-g)*Pstar[i]*Qstar*r1_Pr2_Q;    
        grad[nfact + 1] += difexpg*Qstar*r1_Pr2_Q;
        grad[nfact + 2] += difexpu*Pstar[i]*r1_Pr2_Q;
        for(int j = 0; j < nfact; ++j)
            grad[j] += Theta(i, j)*Pstar[i]*Qstar*(u-g)*r1_Pr2_Q;
        if(estHess){
            double Pstar2 = Pstar[i]*Pstar[i]; 
            double Pstar3 = Pstar[i]*Pstar[i]*Pstar[i];
            hess(nfact,nfact) = hess(nfact,nfact) + (r1_P * (ugD * (Pstar[i] - 3*Pstar2 + 2*Pstar3)) -
                                              r1_P2 * (ugD * (Pstar[i] - Pstar2))*(ugD * (Pstar[i] - Pstar2)) +
                                              r2_Q * (ugD * (-Pstar[i] + 3*Pstar2 - 2*Pstar3)) -
                                              r2_Q2 * (ugD * (-Pstar[i] + Pstar2))*(ugD * (-Pstar[i] + Pstar2)));
            hess(gloc,gloc) = hess(gloc,gloc) + r1_P * (g_1g * (2.0*gm1 - 1.0 - 2.0*gm1*Pstar[i] + Pstar[i])) - 
                                  r1_P2 * (g_1g * (1.0 - Pstar[i])) * (g_1g * (1.0 - Pstar[i])) +
                                  r2_Q * (g_1g * (-2.0*gm1 + 1.0 + 2.0*gm1*Pstar[i] - Pstar[i])) - 
                                  r2_Q2 * (g_1g * (-1.0 + Pstar[i])) * (g_1g * (-1.0 + Pstar[i]));
            hess(uloc,uloc) = hess(uloc,uloc) + r1_P * (2.0*u_1u*um1*Pstar[i]) - r1_P * (u_1u*Pstar[i]) - r1_P2 *(u_1u*u_1u*Pstar2) - 
                                r2_Q * (2.0*u_1u*um1*Pstar[i]) + r2_Q * (u_1u*Pstar[i]) - r2_Q2 *(u_1u*u_1u*Pstar2);
            hess(nfact, gloc) = hess(nfact, gloc)+ (r1_P * (g_1g * (-Pstar[i] + Pstar2)) -
                         r1_P2 * (ugD * (Pstar[i] - Pstar2)) * g_1g * Qstar +
                         r2_Q * (g_1g * (Pstar[i] - Pstar2)) -
                         r2_Q2 * (ugD * (-Pstar[i] + Pstar2)) * g_1g * -Qstar);
            hess(nfact, uloc) = hess(nfact, uloc) + (r1_P * (u_1u * (Pstar[i] - Pstar2)) -
                         r1_P2 * (ugD * (Pstar[i] - Pstar2)) * u_1u * Pstar[i] +
                         r2_Q * (u_1u * (-Pstar[i] + Pstar2)) +
                         r2_Q2 * (ugD * (-Pstar[i] + Pstar2)) * u_1u * Pstar[i]);
            hess(gloc, uloc) = hess(gloc, uloc) + -r1_P2 * (g_1g * (1.0 - Pstar[i])) * u_1u * Pstar[i] + 
                                    r2_Q2 * (g_1g * (-1.0 + Pstar[i])) * u_1u * Pstar[i];
            for(int k = 0; k < nfact; ++k){
                for(int j = 0; j < nfact; ++j){
                    if(k <= j){
                        hess(k, j) = hess(k, j) + (r1_P * (ugD * Theta(i,k) * Theta(i,j) * (Pstar[i] - 3*Pstar2 + 2*Pstar3)) -
                                               r1_P2 * (ugD * Theta(i,k) * (Pstar[i] - Pstar2)) *
                                                  (ugD * Theta(i,j) * (Pstar[i] - Pstar2)) +
                                               r2_Q * (ugD * Theta(i,k) * Theta(i,j) *
                                                   (-Pstar[i] + 3*Pstar2 - 2*Pstar3)) -
                                               r2_Q2 * (ugD * Theta(i,k) * (-Pstar[i] + Pstar2)) *
                                                  (ugD * Theta(i,j) * (-Pstar[i] + Pstar2)));
                    }
                }
            }
            for(int k = 0; k < nfact; ++k){
                hess(k, nfact) = hess(k, nfact) +
                    (r1_P * (ugD * Theta(i,k) * (Pstar[i] - 3*Pstar2 + 2*Pstar3)) -
                             r1_P2 * (ugD * Theta(i,k) * (Pstar[i] - Pstar2)) *
                                (ugD * (Pstar[i] - Pstar2)) +
                             r2_Q * (ugD * Theta(i,k) * (-Pstar[i] + 3*Pstar2 - 2*Pstar3)) -
                             r2_Q2 * (ugD * Theta(i,k) * (-Pstar[i] + Pstar2)) *
                             (ugD * (-Pstar[i] + Pstar2)));
                hess(k, gloc) = hess(k, gloc) +
                    (r1_P * (g_1g * Theta(i,k) * (-Pstar[i] + Pstar2)) -
                             r1_P2 * (ugD * Theta(i,k) * (Pstar[i] - Pstar2)) * g_1g * Qstar +
                             r2_Q * (g_1g * Theta(i,k) * (Pstar[i] - Pstar2)) -
                             r2_Q2 * (ugD * Theta(i,k) * (-Pstar[i] + Pstar2) ) * g_1g * (Pstar[i] - 1.0));
                hess(k, uloc) = hess(k, uloc) +  
                    (r1_P * (u_1u * Theta(i,k) * (Pstar[i] - Pstar2)) -
                             r1_P2 * (ugD * Theta(i,k) * (Pstar[i] - Pstar2)) * u_1u * Pstar[i] +
                             r2_Q * (u_1u * Theta(i,k) * (-Pstar[i] + Pstar2)) +
                             r2_Q2 * (ugD * Theta(i,k) * (-Pstar[i] + Pstar2) ) * u_1u * Pstar[i]);
            }
        }
    }
    if(estHess){
        for(int i = 0; i < hess.nrow(); ++i)
            for(int j = 0; j < hess.ncol(); ++j)
                if(i > j)
                    hess(i, j) = hess(j,i);
    }
    
    List ret;    
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
	return(ret);
	END_RCPP
}

