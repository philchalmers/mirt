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
    int ncol = dat.ncol();
    for(int CAT = 0; CAT < ncol; CAT++){
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
    int ncol = dat.ncol();
    for(int CAT = 0; CAT < ncol; CAT++){
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
        SEXP Rnumak2D2, SEXP RnumakDTheta_numsum) 
{		
    BEGIN_RCPP
    IntegerVector Pnfact(Rnfact), Pncat(Rncat), Pakind(Rakind), Pdind(Rdind);
    int nfact = Pnfact(0), ncat = Pncat(0), akind = Pakind(0), dind = Pdind(0) - 2,
        i, j, k, n, N;
    NumericVector a(Ra), ak(Rak), d(Rd), PD(RD), Prior(RPrior), ak2(Rak2),
                  numsum(Rnumsum), numakD(RnumakD), numak2D2(Rnumak2D2), 
                  aTheta(RaTheta), aTheta2(RaTheta2), dL(nfact + ncat*2);
    NumericVector unitNvec(aTheta.length()); 
    unitNvec.fill(1.0);
    double D = PD(0), D2 = PD(0)*PD(0), tmp;
    NumericMatrix Theta(RTheta), P(RP), num(Rnum), P2(RP2), P3(RP3),
                  dat_num(Rdat_num), numakDTheta_numsum(RnumakDTheta_numsum), 
                  d2L(nfact + ncat*2, nfact + ncat*2), dat(Rdat);
    N = dat.nrow();
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
        dL(dind - i) = sumPrior(tmpvec2, Prior);
    }

    //hess
    //a's
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
                d2L(j, dind - k) = sumPrior(tmpvec2, Prior);
                d2L(dind - k, j) = d2L(j, dind - k);
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
        d2L(dind - j, dind - j) = sumPrior(tmpvec2, Prior);
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
                d2L(dind - i, dind - j) = sumPrior(tmpvec2, Prior); 
                d2L(dind - j, dind - i) = d2L(dind - i, dind - j);
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
                d2L(dind - j, akind + i) = sumPrior(tmpvec, Prior);
                d2L(akind + i, dind - j) = d2L(dind - j, akind + i);
            } else {
                offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta, D2, i);
                for(n = 0; n < N; n++){
                    tmpvec(n) = dat_num(n,i) * (-D2*aTheta(n)*P(n,i)*P(n,j) + 2*P2(n,i) *D2*aTheta(n)*P(n,j)) * numsum(n) + 
                        dat_num(n,i) * (D*P(n,i) - P2(n,i) * D) * D * aTheta(n) * num(n,j) + offterm(n);                
                }
                d2L(akind + i, dind - j) = sumPrior(tmpvec, Prior);
                d2L(dind - j, akind + i) = d2L(akind + i, dind - j);
            }
        }
    }            
    
 

    List ret;
    ret["grad"] = dL;
    ret["hess"] = d2L;
	return(ret);
	END_RCPP
}



/*
  .Call('dparsNominal', a, ak, d, Theta, D, Prior, P, num, dat, nfact, ncat, akind,
  dind, ak2, P2, P3, aTheta, dat_num, numsum, numakD, numakDTheta_numsum)
 */
