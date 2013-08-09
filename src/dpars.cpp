#include"Misc.h"

static NumericVector makeOffterm(const NumericMatrix &dat, const NumericVector &p, const NumericVector &aTheta, 
        const int &cat)
{
    NumericVector ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); CAT++){
        if(CAT == cat) continue;
        for(int n = 0; n < ret.length(); n++)
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
        for(int n = 0; n < ret.length(); n++)
            ret(n) += dat(n, CAT) * p1(n) * p2(n) * aTheta(n);
    }
    return(ret);
}

RcppExport SEXP dparsNominal(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta,
        SEXP RP, SEXP Rnum, SEXP Rdat, SEXP Rnfact, SEXP Rncat,
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
    NumericVector a(Ra), ak(Rak), d(Rd), ak2(Rak2),
                  numsum(Rnumsum), numakD(RnumakD), numak2D2(Rnumak2D2), 
                  aTheta(RaTheta), aTheta2(RaTheta2), dL(nfact + ncat*2);
    NumericVector unitNvec(aTheta.length()); 
    unitNvec.fill(1.0);    
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
                tmpvec(n) += dat_num(n,i)*(ak(i)*Theta(n,j)*P(n,i) - 
                        P(n,i)*numakDTheta_numsum(n,j))*numsum(n);
            }
        }
        dL(j) = sum(tmpvec);               
    }
    for(i = 0; i < ncat; i++){ 
        offterm = makeOffterm(dat, P(_,i), aTheta,  i);
        offterm2 = makeOffterm(dat, P(_,i), unitNvec,  i);
        for(n = 0; n < N; n++){ 
            tmpvec(n) = dat_num(n,i)*(aTheta(n)*P(n,i) - P2(n,i)*aTheta(n))*numsum(n) - offterm(n);
            tmpvec2(n) = dat_num(n,i)*(P(n,i) - P2(n,i))*numsum(n) - offterm2(n);
        }
        dL(akind + i) = sum(tmpvec);
        dL(dind + i) = sum(tmpvec2);
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
        for(j = 0; j < nfact; j++){
            for(k = 0; k < ncat; k++){
                tmpvec.fill(0.0);
                tmpvec2.fill(0.0);
                for(i = 0; i < ncat; i++){
                    for(n = 0; n < N; n++){
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
        for(j = 0; j < ncat; j++){
            tmpvec = makeOffterm(dat, P2(_,j), aTheta2, j);
            tmpvec2 = makeOffterm(dat, P(_,j), aTheta2, j);
            for(n = 0; n < N; n++)
                offterm(n) = tmpvec(n) - tmpvec2(n);
            tmpvec = makeOffterm(dat, P2(_,j), unitNvec, j);
            tmpvec2 = makeOffterm(dat, P(_,j), unitNvec, j);
            for(n = 0; n < N; n++)
                offterm2(n) = tmpvec(n) - tmpvec2(n);
            for(n = 0; n < N; n++){
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
            for(i = 0; i < ncat; i++){
                if(j < i){   
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta2, i);
                    offterm2 = makeOffterm2(dat, P(_,j), P(_,i), unitNvec, i);
                    for(n = 0; n < N; n++){
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
                    for(n = 0; n < N; n++){
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
                    for(n = 0; n < N; n++){
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

	int i, j, k; 
	NumericMatrix prob(Rprob);
	NumericMatrix Thetas(RThetas);    
    NumericMatrix dat(Rdat);
    IntegerVector Pnzeta(Rnzeta);
    IntegerVector estHess(RestHess);        
    const int nzeta = Pnzeta(0);
    const int nfact = Thetas.ncol();
    const int N = Thetas.nrow();
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
    
    //reorder 
    NumericVector grad(dL.length());
    NumericMatrix hess(d2L.ncol(), d2L.ncol());
    
    for(i = 0; i < nfact; i++)
        grad(i) = dL(i+nzeta);
    for(i = 0; i < nzeta; i++)
        grad(i+nfact) = dL(i);    
    if(estHess(0)){
        for(i = 0; i < nfact; i++){
            for(j = 0; j < nfact; j++){
                hess(i,j) = d2L(i+nzeta, j+nzeta);
            }
        }
        for(i = 0; i < nzeta; i++){
            for(j = 0; j < nzeta; j++){
                hess(i+nfact,j+nfact) = d2L(i, j);
            }
        }
        for(i = 0; i < nfact; i++){
            for(j = 0; j < nzeta; j++){
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

RcppExport SEXP dparsDich(SEXP Rx, SEXP RTheta, SEXP RestHess, SEXP REM, SEXP Rot) 
{		
    BEGIN_RCPP
    
    int i, j;   
    NumericVector P, Pstar, Q, Qstar, r1, r2;	
	
	S4 x(Rx);
	NumericMatrix Theta(RTheta);        
	IntegerVector estHess(RestHess);                
	IntegerVector EM(REM);	
	NumericVector ot(Rot);
	
	NumericVector par = x.slot("par");		
	const int nfact = Theta.ncol();    
	NumericVector a(nfact);
	for(i = 0; i < nfact; i++)
		a(i) = par(i);
	const double d = par(nfact);
    const double g = par(nfact+1);
    const double u = par(nfact+2);    
    const double g0 = 0.0;
    const double u1 = 1.0;     	 
	if(EM(0)){
		NumericMatrix dat = x.slot("rs");
		r1 = dat(_,1);
		r2 = dat(_,0);
	} else {
		NumericMatrix dat = x.slot("dat");
		r1 = dat(_,1);
		r2 = dat(_,0);	
	}	
	
    NumericMatrix hess (nfact + 3, nfact + 3);
    NumericVector grad (nfact + 3);
    NumericVector r1_P, r1_P2, r2_Q2, r2_Q;    
    
    P = itemTrace(a, &d, Theta, &g, &u, ot);
    Pstar = itemTrace(a, &d, Theta, &g0, &u1, ot);
    Q = 1.0 - P;
    Qstar = 1.0 - Pstar;        
    r1_P = r1/P;
    r1_P2 = r1/(P*P);
    r2_Q = r2/Q; 
    r2_Q2 = r2/(Q*Q);
    grad(nfact) = sum((u-g)*Pstar*Qstar*(r1_P - r2_Q));    
    grad(nfact + 1) = sum(Qstar*(r1_P - r2_Q));
    grad(nfact + 2) = sum(Pstar*(r1_P - r2_Q));
    for(i = 0; i < nfact; i++)
        grad(i) = sum(Theta(_, i) * Pstar * Qstar * (u-g) * (r1_P - r2_Q));
        
    if(estHess(0)){
        int gloc = nfact+1; 
        int uloc = nfact+2;
        double ugD2 = (u-g);
        double ugD = (u-g); 
        NumericVector Pstar2 = Pstar*Pstar; 
        NumericVector Pstar3 = Pstar*Pstar*Pstar;
        hess(nfact,nfact) = sum((r1_P * (ugD2 * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                          r1_P2 * (ugD * (Pstar - Pstar2))*(ugD * (Pstar - Pstar2)) +
                                          r2_Q * (ugD2 * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                          r2_Q2 * (ugD * (-Pstar + Pstar2))*(ugD * (-Pstar + Pstar2))));
        hess(gloc,gloc) = -1.0 * sum(Qstar*Qstar *(r1_P2 + r2_Q2));
        hess(uloc,uloc) = -1.0 * sum(Pstar2 *(r1_P2 + r2_Q2));
        hess(nfact, gloc) = sum((r1_P * ((-Pstar + Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Qstar +
                     r2_Q * ((Pstar - Pstar2)) -
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * -Qstar));
        hess(gloc, nfact) = hess(nfact, gloc);
        hess(nfact, uloc) = sum((r1_P * ((Pstar - Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Pstar +
                     r2_Q * ((-Pstar + Pstar2)) +
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * Pstar));
        hess(uloc, nfact) = hess(nfact, uloc);
        hess(gloc, uloc) = sum((-1.0 * r1_P2 * Pstar * Qstar + r2_Q2 * Pstar * (-1.0 + Pstar )));
        hess(uloc, gloc) = hess(gloc, uloc);
        for(i = 0; i < nfact; i++){
            for(j = 0; j < nfact; j++){
                if(i <= j){
                    hess(i, j) = sum((r1_P * (ugD2 * Theta(_,i) * Theta(_,j) * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                           r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) *
                                              (ugD * Theta(_,j) * (Pstar - Pstar2)) +
                                           r2_Q * (ugD2 * Theta(_,i) * Theta(_,j) *
                                               (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                           r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2)) *
                                              (ugD * Theta(_,j) * (-Pstar + Pstar2))));                    
                    hess(j, i) = hess(i, j);
                }
            }
        }
        for(i = 0; i < nfact; i++){
            hess(i, nfact) = 
                sum((r1_P * (ugD2 * Theta(_,i) * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) *
                            (ugD * (Pstar - Pstar2)) +
                         r2_Q * (ugD2 * Theta(_,i) * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2)) *
                         (ugD * (-Pstar + Pstar2))));
            hess(nfact, i) = hess(i, nfact);
            hess(i, gloc) = 
                sum((r1_P * (Theta(_,i) * (-Pstar + Pstar2)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) * Qstar +
                         r2_Q * (Theta(_,i) * (Pstar - Pstar2)) -
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2) ) * (Pstar - 1.0)));
            hess(gloc, i) = hess(i, gloc);
            hess(i, uloc) =  
                sum((r1_P * (Theta(_,i) * (Pstar - Pstar2)) -
                         r1_P2 * (ugD * Theta(_,i) * (Pstar - Pstar2)) * Pstar +
                         r2_Q * (Theta(_,i) * (-Pstar + Pstar2)) +
                         r2_Q2 * (ugD * Theta(_,i) * (-Pstar + Pstar2) ) * Pstar));
            hess(uloc, i) = hess(i, uloc);
        }
    }    
    List ret;    
    ret["grad"] = grad;
    ret["hess"] = hess;
	return(ret);
	END_RCPP
}
