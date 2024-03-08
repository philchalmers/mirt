#include "Misc.h"
#include "traceLinePts.h"
#include "Estep.h"
#include "ggum_derivs.h"

static vector<double> makeOffterm(const NumericMatrix &dat, const NumericVector &p, const vector<double> &aTheta,
        const int &cat)
{
    vector<double> ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); ++CAT){
        if(CAT == cat) continue;
        for(int n = 0; n < dat.nrow(); ++n)
            ret[n] += dat(n, CAT) * p(n) * aTheta[n];
    }
    return(ret);
}

static vector<double> makeOffterm2(const NumericMatrix &dat, const NumericVector &p1, const NumericVector &p2,
        const vector<double> &aTheta, const int &cat)
{
    vector<double> ret(dat.nrow());
    for(int CAT = 0; CAT < dat.ncol(); ++CAT){
        if(CAT == cat) continue;
        for(int n = 0; n < dat.nrow(); ++n)
            ret[n] += dat(n, CAT) * p1(n) * p2(n) * aTheta[n];
    }
    return(ret);
}

static double difexp(const double *x)
{
    double ret = *x * (1.0 - *x);
    return(ret);
}

static void add2outer(NumericMatrix &out, const vector<double> &vec, const double &r)
{
    int sz = vec.size();
    for(int i = 0; i < sz; ++i)
        for(int j = 0; j < sz; ++j)
            out(i,j) = out(i,j) + vec[i] * vec[j] * r;
}

// static void add2hess(NumericMatrix &out, const NumericMatrix &in, const double &r)
// {
//     for(int i = 0; i < in.ncol(); ++i)
//         for(int j = 0; j < in.ncol(); ++j)
//             out(i,j) = out(i,j) + in(i,j) * r;
// }

static void matrixMult(vector<double> &c, const vector<double> &a, const vector<double> &b,
                       const int *dim)
{
    NumericMatrix A(*dim, *dim), B(*dim, *dim), C(*dim, *dim);
    int k = 0;

    for (int j = 0; j < *dim; ++j){
        for (int i = 0; i < *dim; ++i){
            A(i,j) = a[k];
            ++k;
        }
    }
    k = 0;
    for (int j = 0; j < *dim; ++j){
        for (int i = 0; i < *dim; ++i){
            B(i,j) = b[k];
            ++k;
        }
    }
    for (int i = 0; i < *dim; ++i){
        for (int j = 0; j < *dim; ++j) {
            C(i,j) = 0;
            for (k = 0; k < *dim; ++k)
                C(i,j) += A(i,k) * B(k,j);
        }
    }
    k = 0;
    for (int j = 0; j < *dim; ++j) {
        for (int i = 0; i < *dim; ++i){
            c[k] = C(i,j);
            ++k;
        }
    }
}

static void matrixMult4(vector<double> &e, const vector<double> &a, const vector<double> &b,
                        const vector<double> &c, const vector<double> &d, const int *dim)
{
    vector<double> tmp1(*dim * (*dim)), tmp2(*dim * (*dim));
    matrixMult(tmp1, a, b, dim);
    matrixMult(tmp2, tmp1, c, dim);
    matrixMult(e, tmp2, d, dim);
}


static double tr(vector<double> &a, const int *dim)
{
    double trace = 0.0;
    int k = 0;

    for(int j = 0; j < *dim; ++j){
        for(int i = 0; i < *dim; ++i){
            if(i == j)
                trace += a[k];
            ++k;
        }
    }
    return trace;
}

static void matrixSub(vector<double> &c, const vector<double> &a, const vector<double> &b,
                      const int *dim)
{
    for(int i = 0; i < *dim*(*dim); ++i)
        c[i] = a[i] - b[i];
}

static void outer(vector<double> &c, const vector<double> &a, const vector<double> &b,
                  const int *dim)
{
    int k = 0;
    for(int i = 0; i < *dim; ++i){
        for(int j = 0; j < *dim; ++j){
            c[k] = a[j] * b[i];
            ++k;
        }
    }
}

static double inner(vector<double> &a, const vector<double> &b, const vector<double> &c,
                    const int *dim)
{
    int k = 0;
    NumericMatrix B(*dim, *dim);
    double ret = 0.0;
    vector<double> tmp(*dim);

    for(int i = 0; i < *dim; ++i){
        tmp[i] = 0.0;
        for(int j = 0; j < *dim; ++j){
            B(j,i) = b[k];
            ++k;
        }
    }
    for(int i = 0; i < *dim; ++i){
        for(int j = 0; j < *dim; ++j){
            tmp[i] += a[j] * B(j,i);
            ++k;
        }
    }
    for(int i = 0; i < *dim; ++i)
        ret += tmp[i] * c[i];
    return ret;
}

static void symMat(vector<double> &dsig, const int *nfact)
{
    int k = 0;
    NumericMatrix tmp(*nfact, *nfact);

    for(int i = 0; i < *nfact; ++i){
        for(int j = 0; j < *nfact; ++j){
            tmp(i,j) = dsig[k];
            ++k;
        }
    }
    for(int i = 0; i < *nfact; ++i)
        for(int j = 0; j < *nfact; ++j)
            if(i < j)
                tmp(j,i) = tmp(i,j);
    k = 0;
    for(int i = 0; i < *nfact; ++i){
        for(int j = 0; j < *nfact; ++j){
            dsig[k] = tmp(i,j);
            ++k;
        }
    }
}



static void _dgroup(vector<double> &grad, NumericMatrix &hess, const NumericVector &par,
	const NumericMatrix &Theta, const bool &estHess, const bool &randeff)
{
    const int N = Theta.nrow();
    const int nfact = Theta.ncol();
    const int npars = nfact + nfact * (nfact + 1);
    const int nsig = npars - nfact;

    arma::vec mu(nfact);
    arma::mat Sig(nfact, nfact);
    int ind;
    if(randeff){
        ind = 0;
        for (int i = 0; i < nfact; ++i){
            for (int j = 0; j < nfact; ++j){
                if(i <= j){
                    Sig(i,j) = par(ind);
                    Sig(j,i) = Sig(i,j);
                    ++ind;
                }
            }
        }
    } else {
        for (int i = 0; i < nfact; ++i) mu(i) = par(i);
        ind = nfact;
        for (int i = 0; i < nfact; ++i){
            for (int j = 0; j < nfact; ++j){
                if(i <= j){
                    Sig(i,j) = par(ind);
                    Sig(j,i) = Sig(i,j);
                    ++ind;
                }
            }
        }
    }
    //const int npars2 = nfact + nfact * (nfact + 1) / 2;
    arma::mat invSig = inv(Sig);
    arma::vec meanTheta(nfact);
    for (int j = 0; j < nfact; ++j){
        double tmp = 0.0;
        for (int i = 0; i < N; ++i)
            tmp += Theta(i,j) / N;
        meanTheta(j) = tmp;
    }

    arma::mat Dif(N, nfact);
    for (int j = 0; j < nfact; ++j)
        for (int i = 0; i < N; ++i)
            Dif(i,j) = Theta(i,j) - mu(j);
    arma::mat Z = trans(Dif) * Dif;
    arma::vec cMeans = N * (meanTheta - mu);
    arma::vec g1 = invSig * cMeans;
    arma::mat Zdif = Z - (N * Sig);
    arma::mat tmp = invSig * Zdif * invSig;
    for (int j = 0; j < nfact; ++j) tmp(j,j) = tmp(j,j)/2;
    ind = nfact;
    for (int i = 0; i < nfact; ++i){
        grad[i] = g1(i);
        for (int j = 0; j < nfact; ++j){
            if(i <= j){
                grad[ind] = tmp(i,j);
                ++ind;
            }
        }
    }
    if(estHess){
        arma::mat hess2(npars, npars);
        vector<double> invSig2(invSig.begin(), invSig.end());
        const vector<double> cMeans2(cMeans.begin(), cMeans.end());
        const vector<double> Zdif2(Zdif.begin(), Zdif.end());
        vector<double> derv1(npars), derv2(npars), du1(nfact), du2(nfact), dsig1(nsig),
            dsig2(nsig), dZ(nsig), dinvSig2(nsig), tmpmat(nsig), dZdif(nsig), Ndsig2(nsig);
        double s1, s2, s3, s4, s5;

        for(int j = 0; j < npars; ++j){
            for(int i = 0; i < npars; ++i){
                if(i <= j){
                    for(int k = 0; k < npars; ++k){
                        derv1[k] = 0.0;
                        derv2[k] = 0.0;
                    }
                    derv1[i] = 1.0;
                    derv2[j] = 1.0;
                    for(int k = 0; k < nfact; ++k){
                        du1[k] = derv1[k];
                        du2[k] = derv2[k];
                    }
                    for(int k = nfact; k < npars; ++k){
                        dsig1[k-nfact] = derv1[k];
                        dsig2[k-nfact] = derv2[k];
                    }
                    symMat(dsig1, &nfact);
                    symMat(dsig2, &nfact);
                    matrixMult(tmpmat, invSig2, dsig2, &nfact);
                    matrixMult(dinvSig2, tmpmat, invSig2, &nfact);
                    for(int k = 0; k < nsig; ++k)
                        dinvSig2[k] = -1.0 * dinvSig2[k];
                    outer(dZ, cMeans2, du2, &nfact);
                    for(int k = 0; k < nsig; ++k)
                        Ndsig2[k] = N * dsig2[k];
                    matrixSub(dZdif, dZ, Ndsig2, &nfact);
                    matrixMult4(tmpmat, dsig1, dinvSig2, Zdif2, invSig2, &nfact);
                    s1 = 0.5 * tr(tmpmat, &nfact);
                    matrixMult4(tmpmat, dsig1, invSig2, Zdif2, dinvSig2, &nfact);
                    s2 = 0.5 * tr(tmpmat, &nfact);
                    matrixMult4(tmpmat, dsig1, invSig2, dZdif, invSig2, &nfact);
                    s3 = 0.5 * tr(tmpmat, &nfact);
                    s4 = inner(du1, dinvSig2, cMeans2, &nfact);
                    s5 = N * inner(du1, invSig2, du2, &nfact);
                    hess2(i,j) = s1 + s2 + s3 + s4 - s5;
                    hess2(j,i) = hess2(i,j);
                }
            }
        }
        arma::uvec pick(nfact + nfact * (nfact + 1) / 2);
        ind = 0;
        int whichrow = nfact;
        for(int  i = 0; i < nfact; ++i){
            pick(i) = ind;
            ++ind;
        }
        for(int  i = 0; i < nfact; ++i){
            for(int  j = 0; j < nfact; ++j){
                if(i <= j){
                    pick(ind) = whichrow;
                    ++ind;
                }
                ++whichrow;
            }
        }
        arma::mat newh = hess2.submat(pick, pick);
        const int n_cols = newh.n_cols;
        for (int i = 0; i < n_cols; ++i)
            for (int j = 0; j < n_cols; ++j)
                hess(i,j) = newh(i,j);
    }
}

static void _dEta(NumericMatrix &dEta, NumericMatrix &d2Eta, const NumericVector &par,
    const NumericMatrix &Theta, const bool &estHess)
{
    const int nquad = Theta.nrow();
    const int nfact = Theta.ncol();
    const int npars = nfact + nfact * (nfact + 1) / 2;
    NumericMatrix theta(1, nfact);
    vector<double> deta(npars);
    NumericMatrix deta2(npars, npars);

    for(int i = 0; i < nquad; ++i){
        for(int j = 0; j < nfact; ++j)
            theta(0,j) = Theta(i,j);
        _dgroup(deta, deta2, par, theta, estHess, false);
        for(int j = 0; j < npars; ++j)
            dEta(i,j) = deta[j];
        int l = 0;
        for(int j = 0; j < npars; ++j){
            for(int k = j; k < npars; ++k){
                d2Eta(i,l) = deta2(j,k);
                ++l;
            }
        }
    }
}

static void _dgroupEM(vector<double> &grad, NumericMatrix &hess, S4 &obj,
	const NumericMatrix &Theta, const NumericMatrix &itemtrace, const vector<double> &prior,
	const bool &estHess)
{
    NumericVector est = obj.slot("est");
    NumericVector par = obj.slot("par");
    bool ret = true;
    for(int i = 0; i < est.length(); ++i)
        if(est(i)) ret = false;
    if(ret) return;

    const int nquad = Theta.nrow();
    const int nfact = Theta.ncol();
    const int npars2 = nfact + nfact * (nfact + 1) / 2;
    NumericMatrix tabdata = obj.slot("dat");
    const int N = tabdata.nrow();
    const int nitems = tabdata.ncol();

    vector<double> hessvec(npars2*(npars2-1)/2);
    NumericMatrix dEta(nquad, npars2);
    NumericMatrix d2Eta(nquad, npars2*npars2);
    _dEta(dEta, d2Eta, par, Theta, estHess);
    for(int pat = 0; pat < N; ++pat){

        vector<double> L(nquad);
        for(int j = 0; j < nquad; ++j)
            L[j] = prior[j];
        for(int j = 0; j < nquad; ++j){
            for(int i = 0; i < nitems; ++i)
                if(tabdata(pat, i))
                    L[j] *= itemtrace(j, i);
        }
        double denom = 0.0;
        const double maxL = *std::max_element(L.begin(), L.end());
        for(int j = 0; j < nquad; ++j) denom += L[j]/maxL;
        denom *= maxL;
        for(int j = 0; j < nquad; ++j)
            L[j] = L[j]/denom;

        for(int j = 0; j < npars2; ++j){
            double tmp = 0.0;
            for(int k = 0; k < nquad; ++k)
                tmp += (L[k] * dEta(k, j));
            grad[j] += tmp;
        }
        if(estHess){
            for(int j = 0; j < npars2*(npars2-1)/2; ++j){
                double tmp = 0.0;
                for(int k = 0; k < nquad; ++k)
                    tmp += (L[k] * d2Eta(k, j));
                hessvec[j] += tmp;
            }
        }
    }

    if(estHess){
        int k = 0;
        for(int i = 0; i < npars2; ++i){
            for(int j = i; j < npars2; ++j){
                hess(i,j) = hessvec[k];
                hess(j,i) = hess(i,j);
                ++k;
            }
        }
    }
}

static void _dgroupEMCD(vector<double> &grad, NumericMatrix &hess, S4 &obj,
	const NumericMatrix &Theta, const bool &estHess)
{
    NumericVector est = obj.slot("est");
    NumericVector par = obj.slot("par");
    bool ret = true;
    for(int i = 0; i < est.length(); ++i)
        if(est(i)) ret = false;
    if(ret) return;

    const int nquad = Theta.nrow();
    const int nfact = Theta.ncol();
    const int npars2 = nfact + nfact * (nfact + 1) / 2;
    const bool BFACTOR = as<bool>(obj.slot("BFACTOR"));

    if(BFACTOR){
        NumericMatrix rrs = obj.slot("rrs");
        NumericVector rrb = obj.slot("rrb");
        NumericMatrix Thetas = obj.slot("theta");
        NumericMatrix Thetab = obj.slot("Thetabetween");
        const int nsfact = rrs.ncol();
        const int npfact = nfact - nsfact;
        const int nppars = npfact + npfact * (npfact + 1) / 2;
        const int npquad = Thetab.nrow();
        const int nsquad = Thetas.nrow();

        NumericVector parp(nppars);
        IntegerVector bindex = obj.slot("bindex");
        for (int i = 0; i < nppars; ++i)
        	parp(i) = par(bindex(i));
        NumericMatrix dEta(npquad, nppars);
        NumericMatrix d2Eta(npquad, nppars*nppars);
        _dEta(dEta, d2Eta, parp, Thetab, estHess);
        for(int j = 0; j < nppars; ++j){
            double tmp = 0.0;
            for(int k = 0; k < npquad; ++k)
                tmp += (rrb(k) * dEta(k, j));
            grad[bindex(j)] += tmp;
        }
        if(estHess){
            vector<double> hessvec(nppars*(nppars-1)/2);
            for(int j = 0; j < nppars*nppars; ++j){
                double tmp = 0.0;
                for(int k = 0; k < npquad; ++k)
                    tmp += (rrb(k) * d2Eta(k, j));
                hessvec[j] += tmp;
            }
            int k = 0;
            for(int i = 0; i < nppars; ++i){
                for(int j = i; j < nppars; ++j){
                    hess(bindex(i),bindex(j)) = hessvec[k];
                    hess(bindex(j),bindex(i)) = hess(bindex(i),bindex(j));
                    ++k;
                }
            }
        }

        NumericMatrix dEtas(nsquad, 2);
        NumericMatrix d2Etas(nsquad, 3);
        IntegerMatrix sindex = obj.slot("sindex");
        for (int s = 0; s < nsfact; ++s){
        	NumericVector pars(2);
        	for (int i = 0; i < 2; ++i)
        		pars(i) = par(sindex(s, i));
        	_dEta(dEtas, d2Etas, pars, Thetas, estHess);
            for(int j = 0; j < 2; ++j){
                double tmp = 0.0;
                for(int k = 0; k < nsquad; ++k)
                    tmp += (rrs(k, s) * dEtas(k, j));
                grad[sindex(s, j)] += tmp;
            }
            if(estHess){
                vector<double> hessvec(3);
                for(int j = 0; j < 3; ++j){
                    double tmp = 0.0;
                    for(int k = 0; k < nsquad; ++k)
                        tmp += (rrs(k, s) * d2Etas(k, j));
                    hessvec[j] += tmp;
                }
                int k = 0;
                for(int i = 0; i < 2; ++i){
                    for(int j = i; j < 2; ++j){
                        hess(sindex(s, i),sindex(s, j)) = hessvec[k];
                        hess(sindex(s, j),sindex(s, i)) = hess(sindex(s, i),sindex(s, j));
                        ++k;
                    }
                }
            }
        }
    } else {
        NumericVector CD = obj.slot("rr");
        NumericMatrix dEta(nquad, npars2);
        NumericMatrix d2Eta(nquad, npars2*npars2);
        NumericMatrix deta2(npars2, npars2);
        _dEta(dEta, d2Eta, par, Theta, estHess);
        for(int j = 0; j < npars2; ++j){
            double tmp = 0.0;
            for(int k = 0; k < nquad; ++k)
                tmp += (CD(k) * dEta(k, j));
            grad[j] += tmp;
        }
        if(estHess){
            vector<double> hessvec(npars2*npars2);
            for(int j = 0; j < npars2*(npars2+1)/2; ++j){
                double tmp = 0.0;
                for(int k = 0; k < nquad; ++k)
                    tmp += (CD(k) * d2Eta(k, j));
                hessvec[j] += tmp;
            }
            int k = 0;
            for(int i = 0; i < npars2; ++i){
                for(int j = i; j < npars2; ++j){
                    hess(i,j) = hessvec[k];
                    hess(j,i) = hess(i,j);
                    ++k;
                }
            }
        }
    }
}

static void _dgroupLCA(vector<double> &grad, NumericMatrix &hess, S4 &obj,
    const NumericMatrix &Theta, const bool &estHess)
{
    NumericVector est = obj.slot("est");
    bool ret = true;
    for(int i = 0; i < est.length(); ++i)
        if(est(i)) ret = false;
    if(ret) return;

    const int nquad = Theta.nrow();
    const int npars = nquad - 1;
    NumericVector CD = obj.slot("rr");
    NumericVector P = obj.slot("density");
    NumericVector P2(P.length());
    NumericVector P3(P.length());
    for(int i = 0; i < nquad; ++i){
        P2(i) = P(i) * P(i);
        P3(i) = P2(i) * P(i);
    }

    for(int p = 0; p < npars; ++p){
        double g = 0.0;
        for(int i = 0; i < nquad; ++i){
            if(p == i){
                g += CD(i) * (P(i) - P2(i)) / P(i);
            } else {
                g -= CD(i) * P(p);
            }
        }
        grad[p] = g;
    }
    if(estHess){
        for(int p = 0; p < npars; ++p){
            for(int q = p; q < npars; ++q){
                double g = 0.0;
                if(p == q){
                    for(int i = 0; i < nquad; ++i){}

                } else {
                    for(int i = 0; i < nquad; ++i){}
                }
                hess(p,q) = g;
                hess(q,p) = g;
            }
        }
        Rprintf("Hessian for LCA hyper-parameters not defined.\n");
    }
}

RcppExport SEXP dgroup(SEXP Robj, SEXP RTheta, SEXP Ritemtrace, SEXP RestHess, SEXP Rrandeff,
        SEXP REM, SEXP REMcomplete)
{
    S4 obj(Robj);
    NumericMatrix Theta(RTheta);
    NumericMatrix itemtrace(Ritemtrace);
    const bool estHess = as<bool>(RestHess);
    const bool randeff = as<bool>(Rrandeff);
    const bool EM = as<bool>(REM);
    const bool EMcomplete = as<bool>(REMcomplete);
    const int nfact = Theta.ncol();
    const int npars2 = nfact + nfact * (nfact + 1) / 2;

    vector<double> grad(npars2);
    int dim = 0; 
    if(estHess) dim = npars2;
    NumericMatrix hess(dim, dim);
    if(EM){
        if(EMcomplete){
            _dgroupEMCD(grad, hess, obj, Theta, estHess);
        } else {
            _dgroupEM(grad, hess, obj, Theta, itemtrace, grad, estHess);
        }
    } else {
        NumericVector par = obj.slot("par");
        _dgroup(grad, hess, par, Theta, estHess, randeff);
    }

    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;

    return(ret);
}

static inline double CDLL(const vector<double> &par, const NumericMatrix &theta,
	const NumericMatrix &dat, const NumericVector &ot, const int &N, const int &nfact,
	const int &ncat, const int &k, const int &itemclass)
{
	vector<double> P(N*ncat);
	P_switch(P, par, theta, ot, N, ncat, nfact, k, itemclass);
	double LL = 0.0;
	for(int j = 0; j < ncat; ++j)
		for(int i = 0; i < N; ++i)
			LL += dat(i,j) * log(P[i + j*N]);
	return(LL);
}

static inline void _central(vector<double> &grad, NumericMatrix &hess,
    const vector<double> &par, const NumericMatrix &theta,
    const NumericMatrix &dat, const NumericVector &ot, const int &N, const int &nfact,
    const int &ncat, const int &k, const int &itemclass, const bool gradient, const double delta)
{
    const int npar = par.size();
    vector<double> parM(npar);
    for(int i = 0; i < npar; ++i)
        parM[i] = par[i];
    if(gradient){
        for(int i = 0; i < npar; ++i){
            parM[i] = par[i] + delta;
            double U = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
            parM[i] = par[i] - 2*delta;
            double L = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
            grad[i] = (U - L) / (2 * delta);
            parM[i] = par[i];
        }
    } else {
        double delta2 = delta*delta;
        double fx = CDLL(par, theta, dat, ot, N, nfact, ncat, k, itemclass);
        for(int i = 0; i < npar; ++i){
            for(int j = i; j < npar; ++j){
                if(i == j){
                    parM[i] = par[i] + 2*delta;
                    double s1 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    parM[i] = par[i] - 2*delta;
                    double s3 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    hess(i, i) = (s1 - 2 * fx + s3) / (4 * delta2);
                } else {
                    parM[i] = par[i] + delta;
                    parM[j] = par[j] + delta;
                    double s1 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    parM[j] = parM[j] - 2*delta;
                    double s2 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    parM[i] = parM[i] - 2*delta;
                    double s4 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    parM[j] = parM[j] + 2*delta;
                    double s3 = CDLL(parM, theta, dat, ot, N, nfact, ncat, k, itemclass);
                    hess(i, j) = (s1 - s2 - s3 + s4) / (4 * delta2);
                    hess(j, i) = hess(i, j);
                }
                parM[i] = par[i];
                parM[j] = par[j];
            }
        }
    }
}

static inline void mat2vec(vector<double> &ret, const NumericMatrix &mat)
{
    const int ncol = mat.ncol();
    const int nrow = mat.nrow();
    int ind = 0;
    for(int j = 0; j < ncol; ++j){
        for(int i = 0; i < nrow; ++i){
            ret[ind] = mat(i,j);
            ++ind;
        }
    }
}

static inline void _richardson(vector<double> &grad, NumericMatrix &hess,
    const vector<double> &par, const NumericMatrix &theta,
    const NumericMatrix &dat, const NumericVector &ot, const int &N, const int &nfact,
    const int &ncat, const int &k, const int &itemclass, const bool gradient)
{
    const int rr = 4;
    const double four = 4.0;
    if(gradient){
        double delta = .0001;
        const int npar = par.size();
        NumericMatrix R0(npar, rr);
        NumericMatrix R1(npar, rr);
        _central(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, true, delta);
        for(int i = 0; i < npar; ++i) R0(i, 0) = grad[i];
        for(int r = 0; r < (rr-1); ++r){
            delta = delta/2;
            _central(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, true, delta);
            for(int i = 0; i < npar; ++i) R1(i, 0) = grad[i];
            for(int j = 0; j < (r + 1); ++j){
                double jj = j;
                double pwr = pow(four, jj + 1.0);
                for(int i = 0; i < npar; ++i)
                    R1(i, j + 1) = (pwr * R1(i, j) - R0(i, j)) / (pwr - 1);
            }
            for(int j = 0; j < r + 1; ++j)
                for(int i = 0; i < npar; ++i)
                    R0(i,j) = R1(i,j);
        }
        for(int i = 0; i < npar; ++i) grad[i] = R1(i, rr - 1);
    } else {
        double delta = .01;
        const int npar = par.size() * par.size();
        const int nrow = hess.nrow();
        vector<double> dvec(npar);
        NumericMatrix R0(npar, rr);
        NumericMatrix R1(npar, rr);
        _central(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, false, delta);
        mat2vec(dvec, hess);
        for(int i = 0; i < npar; ++i) R0(i, 0) = dvec[i];
        for(int r = 0; r < (rr-1); ++r){
            delta = delta/2;
            _central(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, false, delta);
            mat2vec(dvec, hess);
            for(int i = 0; i < npar; ++i) R1(i, 0) = dvec[i];
            for(int j = 0; j < (r + 1); ++j){
                double jj = j;
                double pwr = pow(four, jj + 1.0);
                for(int i = 0; i < npar; ++i)
                    R1(i, j + 1) = (pwr * R1(i, j) - R0(i, j)) / (pwr - 1);
            }
            for(int j = 0; j < r + 1; ++j)
                for(int i = 0; i < npar; ++i)
                    R0(i,j) = R1(i,j);
        }
        int ind = 0;
        for(int j = 0; j < nrow; ++j){
            for(int i = 0; i < nrow; ++i){
                hess(i, j) = R1(ind, rr - 1);
                ++ind;
            }
        }
        for(int j = 0; j < nrow; ++j){
            for(int i = j; i < nrow; ++i){
                if(i != j){
                    hess(i, j) = (hess(i, j) + hess(j, i))/2;
                    hess(j, i) = hess(i,j);
                }
            }
        }
    }
}

static void d_numerical(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
	const NumericMatrix &theta, const NumericVector &ot, const NumericMatrix &dat,
	const int &N, const int &nfact, const int &ncat, const int &k,
    const int &estHess, const int &itemclass)
{
	const int supported[] = {6, 9, 10, 11, 12}; // supported item class #
	bool run = false;
	for(int i = 0; i < 5; ++i) // length of supported
		if(supported[i] == itemclass) run = true;
	if(!run) return;

    _richardson(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, true);
    if(estHess)
        _richardson(grad, hess, par, theta, dat, ot, N, nfact, ncat, k, itemclass, false);
}

static void d_nominal(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &ncat, const int &israting, const int &estHess)
{
    vector<double> p(N*ncat), pnum(N*ncat);
    P_nominal(p, par, Theta, ot, N, nfact, ncat, 0, israting);
    P_nominal(pnum, par, Theta, ot, N, nfact, ncat, 1, israting);
    const NumericMatrix P = vec2mat(p, N, ncat);
    const NumericMatrix num = vec2mat(pnum, N, ncat);
    NumericMatrix P2(N, ncat), dat_num(N, ncat);
    vector<double> numsum(N, 0.0);
    for(int i = 0; i < N; ++i){
        long double tmpnumsum = 0.0;
        for(int j = 0; j < ncat; ++j){
            P2(i,j) = P(i,j) * P(i,j);
            dat_num(i,j) = dat(i,j) / num(i,j);
            tmpnumsum += num(i,j);
        }
        numsum[i] = tmpnumsum;
    }
    const int akind = nfact;
    const int dind = nfact + ncat;
    vector<double> a(nfact), ak(ncat), d(ncat), ak2(ncat);
    for(int i = 0; i < nfact; ++i)
        a[i] = par[i];
    for(int i = 0; i < ncat; ++i){
        ak[i] = par[i + nfact];
        ak2[i] = ak[i] * ak[i];
        if(israting){
            if(i) d[i] = par[i + nfact + ncat] + par[par.size()-1];
        } else {
            d[i] = par[i + nfact + ncat];
        }
    }
    vector<double> unitNvec(N, 1.0), aTheta(N);
    vector<double> numakD(N);
    NumericMatrix numakDTheta_numsum(N, nfact);
    for(int i = 0; i < N; ++i){
        double tmpa = 0.0;
        for(int j = 0; j < ncat; ++j){
            numakD[i] += num(i,j) * ak[j];
        }
        for(int j = 0; j < nfact; ++j){
            tmpa += a[j] * Theta(i, j);
            numakDTheta_numsum(i, j) = numakD[i] * Theta(i, j) / numsum[i];
        }
        aTheta[i] = tmpa;
    }
    vector<double> tmpvec(N), tmpvec2(N), offterm(N), offterm2(N);

    //grad
    for(int j = 0; j < nfact; ++j){
        std::fill(tmpvec.begin(), tmpvec.end(), 0.0);
        for(int i = 0; i < ncat; ++i){
            for(int n = 0; n < N; ++n){
                tmpvec[n] += dat_num(n,i)*P(n,i)*(ak[i]*Theta(n,j) -
                        numakDTheta_numsum(n,j))*numsum[n];
            }
        }
        grad[j] = vecsum(tmpvec);
    }
    for(int i = 0; i < ncat; ++i){
        offterm = makeOffterm(dat, P(_,i), aTheta, i);
        offterm2 = makeOffterm(dat, P(_,i), unitNvec, i);
        for(int n = 0; n < N; ++n){
            tmpvec[n] = dat_num(n,i)*P(n,i)*aTheta[n]*(1.0 - P(n,i))*numsum[n] - offterm[n];
            tmpvec2[n] = dat_num(n,i)*P(n,i)*(1.0 - P(n,i))*numsum[n] - offterm2[n];
        }
        grad[akind + i] = vecsum(tmpvec);
        grad[dind + i] = vecsum(tmpvec2);
    }

    //hess
    //a's
    if(estHess){
        vector<double> numak2D2(N), aTheta2(N);
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < ncat; ++j){
                numak2D2[i] += num(i,j) * ak2[j];
            }
            aTheta2[i] = aTheta[i] * aTheta[i];
        }

        for(int j = 0; j < nfact; ++j){
            for(int k = 0; k < nfact; ++k){
                if(j <= k){
                    std::fill(tmpvec.begin(), tmpvec.end(), 0.0);
                    for(int i = 0; i < ncat; ++i){
                        for(int n = 0; n < N; ++n){
                            tmpvec[n]+= dat_num(n,i)*P(n,i)*(ak2[i]*Theta(n,j)*Theta(n,k) -
                                    ak[i]*Theta(n,j)*numakDTheta_numsum(n,k) -
                                    ak[i]*Theta(n,k)*numakDTheta_numsum(n,j) +
                                    2*numakD[n]*Theta(n,j)*numakD[n]*Theta(n,k)/ (numsum[n]*numsum[n]) -
                                    numak2D2[n]*Theta(n,j)*Theta(n,k)/numsum[n]) * numsum[n] -
                                dat_num(n,i)*P(n,i)*(ak[i]*Theta(n,j) - numakDTheta_numsum(n,j)) *
                                numsum[n]*ak[i]*Theta(n,k) +
                                dat_num(n,i)*P(n,i)*(ak[i]*Theta(n,j) - numakDTheta_numsum(n,j)) *
                                numakD[n]*Theta(n,k);
                        }
                    }
                    hess(j,k) = vecsum(tmpvec);
                    hess(k, j) = hess(j,k);
                }
            }
        }
        //a's with ak and d
        for(int j = 0; j < nfact; ++j){
            for(int k = 0; k < ncat; ++k){
                std::fill(tmpvec.begin(), tmpvec.end(), 0.0);
                std::fill(tmpvec2.begin(), tmpvec2.end(), 0.0);
                for(int i = 0; i < ncat; ++i){
                    for(int n = 0; n < N; ++n){
                        if(i == k){
                            tmpvec[n] += dat_num(n,i)*P(n,i)*(ak[i]*Theta(n,j)*aTheta[n] -
                                        aTheta[n]*numakDTheta_numsum(n,j) +
                                        Theta(n,j) - 2*ak[i]*Theta(n,j)*aTheta[n]*P(n,i) +
                                        2*aTheta[n]*P(n,i)*numakDTheta_numsum(n,j) -
                                        Theta(n,j)*P(n,i))*numsum[n] -
                                dat_num(n,i)*P(n,i)*aTheta[n]*(1.0 - P(n,i))*numsum[n]*ak[i]*Theta(n,j) +
                                dat_num(n,i)*P(n,i)*aTheta[n]*(1.0 - P(n,i))*(numakD[n]*Theta(n,j));
                        tmpvec2[n] += dat_num(n,i)*P(n,i)*(ak[i]*Theta(n,j) -
                                                            2*ak[i]*Theta(n,j)*P(n,i) -
                                                            numakDTheta_numsum(n,j) +
                                                            2*P(n,i)*numakDTheta_numsum(n,j))*numsum[n] -
                            dat_num(n,i)*P(n,i)*(1.0 - P(n,i))*numsum[n]*ak[i]*Theta(n,j) +
                            dat_num(n,i)*P(n,i)*(1.0 - P(n,i))*(numakD[n]*Theta(n,j));
                        } else {
                            tmpvec[n] += dat(n,i)*P(n,k)*(-ak[k]*aTheta[n]*Theta(n,j) +
                                aTheta[n]*numakDTheta_numsum(n,j) - Theta(n,j));
                            tmpvec2[n] += dat(n,i)*P(n,k)*(-ak[k]*Theta(n,j) +
                                numakDTheta_numsum(n,j));
                        }
                    }
                    hess(j, akind + k) = vecsum(tmpvec);
                    hess(akind + k, j) = hess(j, akind + k);
                    hess(j, dind + k) = vecsum(tmpvec2);
                    hess(dind + k, j) = hess(j, dind + k);
                }
            }
        }
        // //ak's and d's
        for(int j = 0; j < ncat; ++j){
            tmpvec = makeOffterm2(dat, P(_,j), P(_,j), aTheta2, j);
            tmpvec2 = makeOffterm(dat, P(_,j), aTheta2, j);
            for(int n = 0; n < N; ++n)
                offterm[n] = tmpvec[n] - tmpvec2[n];
            tmpvec = makeOffterm2(dat, P(_,j), P(_,j), unitNvec, j);
            tmpvec2 = makeOffterm(dat, P(_,j), unitNvec, j);
            for(int n = 0; n < N; ++n)
                offterm2[n] = tmpvec[n] - tmpvec2[n];
            for(int n = 0; n < N; ++n){
                tmpvec[n] = P(n,j)*(dat_num(n,j)*aTheta2[n]*(1.0 - 3*P(n,j) + 2*P(n,j)*P(n,j)) * numsum[n] -
                                                dat_num(n,j)*aTheta2[n]*(1.0 - P(n,j))*numsum[n] +
                                                dat(n,j)*aTheta2[n]*(1.0 - P(n,j))) + offterm[n];
                tmpvec2[n] = P(n,j)*dat(n,j)*(1.0/num(n,j)*(1.0 - 3*P(n,j) + 2*P(n,j)*P(n,j)) * numsum[n] -
                    1.0/num(n,j)*(1.0 - P(n,j))*numsum[n] +
                    (1.0 - P(n,j))) + offterm2[n];
            }
            hess(akind + j, akind + j) = vecsum(tmpvec);
            hess(dind + j, dind + j) = vecsum(tmpvec2);
            for(int i = 0; i < ncat; ++i){
                if(j < i){
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta2, i);
                    offterm2 = makeOffterm2(dat, P(_,j), P(_,i), unitNvec, i);
                    for(int n = 0; n < N; ++n){
                        tmpvec[n] = dat_num(n,i) * (-aTheta2[n]*P(n,i)*P(n,j) + 2*P2(n,i) *aTheta2[n]*P(n,j))*numsum[n] +
                                     dat_num(n,i) * (aTheta[n]*P(n,i) - P2(n,i) * aTheta[n])*aTheta[n]*num(n,j)+offterm[n];
                        tmpvec2[n] = dat_num(n,i) * (-P(n,i)*P(n,j) + 2*P2(n,i) *P(n,j)) * numsum[n] +
                            dat_num(n,i) * (P(n,i) - P2(n,i)) * num(n,j) + offterm2[n];
                    }
                    hess(akind + i, akind + j) = vecsum(tmpvec);
                    hess(akind + j, akind + i) = hess(akind + i, akind + j);
                    hess(dind + i, dind + j) = vecsum(tmpvec2);
                    hess(dind + j, dind + i) = hess(dind + i, dind + j);
                }
                if(abs(j-i) == 0){
                    tmpvec = makeOffterm2(dat, P(_,i), P(_,i), aTheta, i);
                    tmpvec2 = makeOffterm(dat, P(_,i), aTheta, i);
                    for(int n = 0; n < N; ++n){
                        offterm[n] = tmpvec[n] - tmpvec2[n];
                        tmpvec[n] = dat_num(n,i)*P(n,i)*aTheta[n]*(1.0 - 3*P(n,i) +
                                2*P2(n,i))*numsum[n] - dat_num(n,i)*aTheta[n]*P(n,i)*(1.0 -
                                P(n,i))*numsum[n] + dat(n,i)*P(n,i)*(1.0 -
                                P(n,i))*aTheta[n] + offterm[n];
                    }
                    hess(dind + j, akind + i) = vecsum(tmpvec);
                    hess(akind + i, dind + j) = hess(dind + j, akind + i);
                } else {
                    offterm = makeOffterm2(dat, P(_,j), P(_,i), aTheta, i);
                    for(int n = 0; n < N; ++n){
                        tmpvec[n] = dat_num(n,i) * (-aTheta[n]*P(n,i)*P(n,j) + 2*P2(n,i) *aTheta[n]*P(n,j)) * numsum[n] +
                            dat_num(n,i) * P(n,i) * (1.0 - P(n,i)) * aTheta[n] * num(n,j) + offterm[n];
                    }
                    hess(akind + i, dind + j) = vecsum(tmpvec);
                    hess(dind + j, akind + i) = hess(akind + i, dind + j);
                }
            }
        }
    }
}

static void d_nominal2(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &ncat, const int &israting, const int &estHess)
{
    const int dind = nfact + ncat*nfact;

    vector<double> p(N*ncat), pnum(N*ncat), Q(N);
    P_nominal2(p, par, Theta, ot, N, nfact, ncat, 0, 0);
    P_nominal2(pnum, par, Theta, ot, N, nfact, ncat, 1, 0);
    const NumericMatrix P = vec2mat(p, N, ncat);
    const NumericMatrix num = vec2mat(pnum, N, ncat);
    NumericMatrix P2(N, ncat);

    for(int i = 0; i < N; ++i){
        double tmp = 0.0;
        for(int k = 0; k < ncat; ++k){
            tmp += num(i, k);
        }
        Q[i] = 1.0 / tmp;
    }

    vector<double> a(nfact), d(ncat);
    NumericMatrix ak(ncat, nfact);
    for(int j = 0; j < nfact; ++j)
        a[j] = par[j];
    int ind = nfact;
    for(int j = 0; j < nfact; ++j){
        for(int k = 0; k < ncat; ++k){
            ak(k, j) = par[ind];
            ++ind;
        }
    }
    for(int k = 0; k < ncat; ++k)
        d[k] = par[k + nfact + ncat*nfact];

    for(int i = 0; i < N; ++i){
      //long double tmpnumsum = 0.0;
        for(int k = 0; k < ncat; ++k){
            P2(i,k) = P(i,k) * P(i,k);
            //tmpnumsum += num(i,j);
        }
        //numsum[i] = tmpnumsum;
    }

    NumericMatrix NumSum(N,nfact);
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < nfact; ++j){
            double tmp = 0.0;
            for(int k = 0; k < ncat; ++k)
                tmp += ak(k,j) * Theta(i,j) * num(i,k);
            NumSum(i,j) = tmp;
        }
    }

    NumericMatrix akakDThetaTheta(N, nfact);
    vector<double> unitNvec(N, 1.0), tmpvec(N), tmpvec2(N), offterm(N), offterm2(N);

    //grad
    ind = 0;
    for(int j = 0; j < nfact; ++j){
        double tmp = 0.0;
        for(int k = 0; k < ncat; ++k){
	    for(int i = 0; i < N; ++i){
      	      tmp += dat(i,k) * (ak(k,j) * Theta(i,j) - Q[i] * NumSum(i,j));
	    }
        }
        grad[ind] = tmp;
        ++ind;
    }
    ind = nfact*ncat + nfact;
    for(int k = 0; k < ncat; ++k){
        offterm2 = makeOffterm(dat, P(_,k), unitNvec, k);
        double tmp = 0.0;
        for (int i = 0; i < N; ++i){
	    tmp += dat(i,k) * (1.0 - P(i,k)) - offterm2[i];
	}
        grad[ind] = tmp;
        ++ind;
    }


    //hess
    if(estHess){

        //a's
        for(int j = 0; j < nfact; ++j){

  	    std::fill(akakDThetaTheta.begin(), akakDThetaTheta.end(), 0.0);
  	    for(int i = 0; i < N; ++i){
  	        for(int d = 0; d < nfact; ++d){
  	            for(int k=0; k < ncat; ++k){
		        akakDThetaTheta(i, d) += P(i,k) * ak(k,d) * Theta(i, d) * ak(k,j) * Theta(i, j);
	            }
	         }
	     }

             for(int d = 0; d < nfact; ++d){
                 if(j <= d){
                     std::fill(tmpvec.begin(), tmpvec.end(), 0.0);

                     for(int i = 0; i < ncat; ++i){
                         for(int n = 0; n < N; ++n){
			   tmpvec[n]+= dat(n,i)*Q[n]*NumSum(n,j)*Q[n]*NumSum(n,d) - dat(n,i)*akakDThetaTheta(n,d);
                         }
                     }
                     hess(j, d) = vecsum(tmpvec);
                     hess(d, j) = hess(j, d);
                 }
             }
         }

         //a's and d's
         for(int j = 0; j < nfact; ++j){
             for(int k = 0; k < ncat; ++k){
                 std::fill(tmpvec2.begin(), tmpvec2.end(), 0.0);
                 for(int i = 0; i < ncat; ++i){
                     for(int n = 0; n < N; ++n){
		       tmpvec2[n] += dat(n,i)*P(n,k)*(-ak(k,j)*Theta(n,j)+Q[n]*NumSum(n,j));
                     }
                     hess(j, dind + k) = vecsum(tmpvec2);
                     hess(dind + k, j) = hess(j, dind + k);
                 }
             }
         }

         // d's
         for(int j = 0; j < ncat; ++j){
             tmpvec = makeOffterm2(dat, P(_,j), P(_,j), unitNvec, j);
             tmpvec2 = makeOffterm(dat, P(_,j), unitNvec, j);
             for(int n = 0; n < N; ++n)
                 offterm2[n] = tmpvec[n] - tmpvec2[n];
             for(int n = 0; n < N; ++n)
	       tmpvec2[n] = dat(n,j)*(P2(n,j)-P(n,j))+offterm2[n];
             hess(dind + j, dind + j) = vecsum(tmpvec2);

             for(int i = 0; i < ncat; ++i){
                 if(j < i){
                     offterm2 = makeOffterm2(dat, P(_,j), P(_,i), unitNvec, i);
                     for(int n = 0; n < N; ++n){
		       tmpvec2[n] = dat(n,i)*(P(n,i)*P(n,j)) + offterm2[n];
                     }
                     hess(dind + i, dind + j) = vecsum(tmpvec2);
                     hess(dind + j, dind + i) = hess(dind + i, dind + j);
                 }
	     }

         }

    }
}

RcppExport SEXP dparsNominal(SEXP Rx, SEXP RTheta, SEXP Roffterm,
    SEXP Risrating, SEXP RestHess)
{
    BEGIN_RCPP

    S4 x(Rx);
    const NumericMatrix dat = x.slot("dat");
    const vector<double> par = as< vector<double> >(x.slot("par"));
    const int ncat = as<int>(x.slot("ncat"));
    const NumericMatrix Theta(RTheta);
    const NumericVector ot(Roffterm);
    const int israting = as<int>(Risrating);
    const int N = Theta.nrow();
    const int nfact = Theta.ncol();
    const int estHess = as<int>(RestHess);
    const int has_mat = as<int>(x.slot("mat"));
    int size = par.size();
    vector<double> grad(size);
    int dim = 0; 
    if(estHess) dim = size;
    NumericMatrix hess(dim, dim);
    if(has_mat){
        d_nominal2(grad, hess, par, Theta, ot, dat, N, nfact, ncat, israting, estHess);
    } else {
        d_nominal(grad, hess, par, Theta, ot, dat, N, nfact, ncat, israting, estHess);
    }

    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
    return(ret);

    END_RCPP
}

void d_poly(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &nzeta, const int &estHess)
{
    vector<double> Pprob(N * (nzeta + 2));
    P_graded(Pprob, par, Theta, ot, N, nfact, nzeta, 0, 0);
    const NumericMatrix prob = vec2mat(Pprob, N, nzeta + 2);

    vector<double> Pk(N), Pk_1(N), Pk_p1(N), PQ_1(N), PQ(N), PQ_p1(N),
            Pk_1Pk(N), Pk_Pkp1(N), dif1(N), dif1sq(N), dif2(N),
            dif2sq(N), tmp1(N), tmp2(N), tmp3(N), csums(nfact);
    NumericMatrix P(N,nzeta+2), PQfull(N,nzeta+2), mattmp(N,nfact), d2Louter;
    NumericVector dL(grad.size());
    NumericMatrix d2L(grad.size(), grad.size());
    vector<int> factind(nfact);
    for(int j = 0; j < (nzeta + 2); ++j){
        for(int i = 0; i < N; ++i){
            P(i,j) = prob(i,j);
            PQfull(i,j) = prob(i,j) * (1.0 - prob(i,j));
        }
    }
    for(int j = 0; j < nfact; ++j)
        factind[j] = nzeta + j;
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
                        mattmp(i,k) = tmp1[i] * Theta(i,k) + tmp2[i] * Theta(i,k) -
                            tmp3[i] * Theta(i,k);
                        csums[k] += mattmp(i,k);
                    }
                }
                for(int i = 0; i < nfact; ++i){
                    d2L(j,factind[i]) = csums[i];
                    d2L(factind[i],j) = csums[i];
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
                mattmp(i,k) = dif1[i] * (PQ_1[i] - PQ[i]) * Theta(i,k);
                csums[k] += mattmp(i,k);
            }
        }
        for(int i = 0; i < nfact; ++i)
            dL(factind[i]) += csums[i];

        if(estHess){
            d2Louter = polyOuter(Theta, Pk, Pk_1, PQ_1, PQ, dif1sq, dif1);
            for(int k = 0; k < nfact; ++k)
                for(int i = 0; i < nfact; ++i)
                    d2L(factind[i],factind[k]) += d2Louter(i,k);
        }
    }

    //reorder
    for(int i = 0; i < nfact; ++i)
        grad[i] = dL(i+nzeta);
    for(int i = 0; i < nzeta; ++i)
        grad[i+nfact] = dL(i);
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
}

RcppExport SEXP dparsPoly(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Rdat, SEXP Rnzeta, SEXP RestHess)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericVector ot(Rot);
	const NumericMatrix Theta(RTheta);
    const NumericMatrix dat(Rdat);
    const int nzeta = as<int>(Rnzeta);
    const int estHess = as<int>(RestHess);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int dim = 0; 
    if(estHess) dim = nfact + nzeta;
    NumericMatrix hess(dim, dim);
    vector<double> grad(nfact + nzeta);
    d_poly(grad, hess, par, Theta, ot, dat, N, nfact, nzeta, estHess);
    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
	return(ret);

	END_RCPP
}

void d_gpcmIRT(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &nzeta, const int &estHess)
{
    if(estHess){
         d_numerical(grad, hess, par, Theta, ot,
            dat, N, nfact, nzeta + 1, 0, estHess, 6);
    }
    vector<double> Pprob(N * (nzeta + 1));
    P_gpcmIRT(Pprob, par, Theta, ot, N, 1, nzeta);
    const NumericMatrix P = vec2mat(Pprob, N, nzeta + 1);
    const int parsize = par.size();
    const int ncat = parsize - 1;

    const double a = par[0];
    vector<double> b(ncat-1), bsum(ncat, 0.0);
    for(int i = 1; i < (parsize-1); ++i){
        b[i-1] = par[i];
        bsum[i] = b[i-1] + bsum[i-1];
    }

    for (int i = 0; i < N; ++i){
        vector<double> r1_P(ncat);
        double psia = 0.0, psic = 0.0;
        for (int j = 0; j < ncat; ++j){
            r1_P[j] = dat(i,j) / P(i,j);
            psia += (j * Theta(i, 0) - bsum[j]) * P(i,j);
            psic += j * P(i,j);
        }

        grad[0] += dat(i,0) * (-psia);
        grad[parsize-1] += dat(i,0) * (-psic);
        for (int j = 1; j < ncat; ++j){
            grad[0] += r1_P[j] * ( (j * Theta(i, 0) - bsum[j]) * P(i,j) - P(i,j) * psia );
            grad[parsize-1] += r1_P[j] * ( j * P(i,j) - P(i,j) * psic );
        }
        for (int j = 0; j < ncat - 1; ++j){
            double psib = 0.0;
            for (int k = j+1; k < ncat; ++k)
                psib += a*P(i,k);
            for (int k = 0; k < j+1; ++k)
                grad[j+1] += dat(i,k) * psib;
            for (int k = j+1; k < ncat; ++k)
                grad[j+1] += r1_P[k] * (-a * P(i,k) + P(i,k) * psib );
        }
    }
}

RcppExport SEXP dparsgpcmIRT(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Rdat, SEXP Rnzeta, SEXP RestHess)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericVector ot(Rot);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix dat(Rdat);
    const int nzeta = as<int>(Rnzeta);
    const int estHess = as<int>(RestHess);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int dim = 0;
    if(estHess) dim = nfact + nzeta;
    NumericMatrix hess(dim, dim);
    vector<double> grad(nfact + nzeta);
    d_gpcmIRT(grad, hess, par, Theta, ot, dat, N, nfact, nzeta, estHess);
    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
    return(ret);

    END_RCPP
}

void d_lca(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericMatrix &item_Q,
    const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &estHess)
{
    const int ncat = dat.ncol();
    if(estHess){
        d_numerical(grad, hess, par, Theta,
            ot, dat, N, nfact, ncat, 0, estHess, 10);
    }
    vector<double> p(N*ncat);
    P_lca(p, par, Theta, item_Q, N, ncat, nfact, 0);
    const NumericMatrix P = vec2mat(p, N, ncat);

    for (int i = 0; i < N; ++i){
        int ind = 0;
        for (int k = 1; k < ncat; ++k){
            for (int j = 0; j < nfact; ++j){
                double val = dat(i, k) * ( P(i, k) * (1.0-P(i, k)) ) / P(i, k);
                for (int kk = 0; kk < ncat; ++kk)
                    if (kk != k)
                        val -= dat(i, kk) * P(i, k);
                val *= Theta(i, j);
                grad[ind] += val * item_Q(k, j);
                ind++;
            }
        }
    }
}

RcppExport SEXP dparslca(SEXP Rx, SEXP RTheta, SEXP Ritem_Q, SEXP RestHess, SEXP Rdat, SEXP Rot)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rx);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix item_Q(Ritem_Q);
    const NumericMatrix dat(Rdat);
    const NumericVector ot(Rot);
    const int estHess = as<int>(RestHess);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int dim = 0;
    if(estHess) dim = par.size();
    NumericMatrix hess (dim, dim);
    vector<double> grad (par.size());
    d_lca(grad, hess, par, Theta, item_Q, ot, dat, N, nfact, estHess);
    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
    return(ret);

    END_RCPP
}

void d_dich(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &estHess)
{
    vector<double> a(nfact);
    for(int i = 0; i < nfact; ++i) a[i] = par[i];
    const double d = par[nfact];
    const double expg = par[nfact+1];
    const double expu = par[nfact+2];
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

    vector<double> P(N), Pstar(N);
    itemTrace(P, Pstar, a, &d, Theta, &g, &u, ot);

    for(int i = 0; i < N; ++i){
        double Q = 1.0 - P[i];
        double Qstar = 1.0 - Pstar[i];
        double r1_P = dat(i, 1) / P[i];
        double r2_Q = dat(i, 0) / Q;
        double r1_Pr2_Q = r1_P - r2_Q;
        grad[nfact] += (u-g)*Pstar[i]*Qstar*r1_Pr2_Q;
        grad[nfact + 1] += difexpg*Qstar*r1_Pr2_Q;
        grad[nfact + 2] += difexpu*Pstar[i]*r1_Pr2_Q;
        for(int j = 0; j < nfact; ++j)
            grad[j] += Theta(i, j)*Pstar[i]*Qstar*(u-g)*r1_Pr2_Q;
        if(estHess){
            double r1_P2 = dat(i, 1) / (P[i]*P[i]);
            double r2_Q2 = dat(i, 0) / (Q*Q);
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
}

RcppExport SEXP dparsDich(SEXP Rx, SEXP RTheta, SEXP RestHess, SEXP Rdat, SEXP Rot)
{
    BEGIN_RCPP

	const vector<double> par = as< vector<double> >(Rx);
    const NumericMatrix Theta(RTheta);
	const NumericMatrix dat(Rdat);
	const NumericVector ot(Rot);
    const int estHess = as<int>(RestHess);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    NumericMatrix hess (nfact + 3, nfact + 3);
    vector<double> grad (nfact + 3);
    d_dich(grad, hess, par, Theta, ot, dat, N, nfact, estHess);
    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
	return(ret);

	END_RCPP
}

void d_ggum(vector<double> &grad, NumericMatrix &hess, const vector<double> &par,
    const NumericMatrix &Theta, const NumericMatrix &dat,
    const int &N, const int &nfact, const int &ncat, const int &estHess)
{
    const int D = nfact;
    const int C = ncat - 1;
    arma::colvec par2(par);
    arma::mat Theta2 = Rcpp::as<arma::mat>(Theta);
    arma::mat Z = Rcpp::as<arma::mat>(dat);

    NumericVector grad_tmp = grad_ggum (par2, Theta2, D, C, Z);
    for (int i = 0; i < grad_tmp.length(); ++i)
        grad[i] = grad_tmp(i);

    if(estHess){
        arma::mat hess_tmp = hess_ggum (par2, Theta2, D, C, Z);
        for(int i = 0; i < hess.nrow(); ++i)
            for(int j = 0; j < hess.ncol(); ++j)
                    hess(i, j) = hess_tmp(j,i);
    }
}

// RcppExport SEXP dparsGGUM(SEXP Rx, SEXP RTheta, SEXP RestHess, SEXP Rdat)
// {
//     BEGIN_RCPP


//     const vector<double> par = as< vector<double> >(Rx);
//     const NumericMatrix Theta(RTheta);
//     const NumericMatrix dat(Rdat);
//     const int estHess = as<int>(RestHess);
//     const int nfact = Theta.ncol();
//     const int N = Theta.nrow();
//     NumericMatrix hess (nfact + 3, nfact + 3);
//     vector<double> grad (nfact + 3);
//     d_ggum(tmpgrad, tmphess, par, theta, dat, N, nfact, ncat, estHess);
//     List ret;
//     ret["grad"] = wrap(grad);
//     ret["hess"] = hess;
//     return(ret);

//     END_RCPP
// }

static void d_priors(vector<double> &grad, NumericMatrix &hess, const int &ind,
    const int &prior_type, const double &prior_1, const double &prior_2, const double &par)
{
    double g = 0.0, h = 0.0;
    if(prior_type == 1){
        const double p2 = prior_2*prior_2;
        g = -(par - prior_1)/p2;
        h = -1.0/p2;
    } else if(prior_type == 2){
        double val = par;
        if(val < 1e-10) val = 1e-10;
        const double lval = log(val);
        const double p2 = prior_2*prior_2;
        const double v2 = val*val;
        g = -(lval - prior_1)/(val * p2) - 1.0/val;
        h = 1.0/v2 - 1.0/(v2 * p2) - (lval - prior_1)/(v2 * p2);
    } else if((prior_type == 3) | (prior_type == 4)){
        double val = par;
        if(prior_type == 4) val = 1 / (1 + exp(-val));
        if(val < 1e-10) val = 1e-10;
        else if(val > 1.0 - 1e-10) val = 1.0 - 1e-10;
        g = (prior_1 - 1.0)/val - (prior_2 - 1.0)/(1.0 - val);
        h = -(prior_1 - 1.0)/(val*val) - (prior_2 - 1.0) / ((1.0 - val) * (1.0 - val));
    }
    grad[ind] += g;
    hess(ind, ind) = hess(ind, ind) + h;
}

static void _computeDpars(vector<double> &grad, NumericMatrix &hess, const List &pars,
    const NumericMatrix &Theta, const NumericMatrix &offterm, const NumericMatrix &itemtrace,
    const vector<double> &prior, const int &nitems, const int &npars,
    const int &estHess, const int &USEFIXED, const int &EM, const bool &EMcomplete,
    const bool &useIprior)
{
    int nfact = Theta.ncol();
    int N = Theta.nrow();
    int has_mat = 0;
    for(int i = 0; i < nitems + EM; ++i){
        S4 item = pars[i];
        int nfact2 = nfact;
        NumericMatrix theta = Theta;
        if(USEFIXED){
            NumericMatrix itemFD = item.slot("fixed.design");
            nfact2 = nfact + itemFD.ncol();
            NumericMatrix NewTheta(Theta.nrow(), nfact2);
            for(int j = 0; j < itemFD.ncol(); ++j)
                NewTheta(_,j) = itemFD(_,j);
            for(int j = 0; j < nfact; ++j)
                NewTheta(_,j + itemFD.ncol()) = Theta(_,j);
            theta = NewTheta;
        }
        vector<double> par = as< vector<double> >(item.slot("par"));
		int par_size = par.size();
        vector<double> tmpgrad(par_size);
        NumericMatrix tmphess(par_size, par_size);
        int itemclass = as<int>(item.slot("itemclass"));
        int ncat;
        if(itemclass > 0)
            ncat = as<int>(item.slot("ncat"));
        int k = 0;
        if(itemclass == 12) k = as<int>(item.slot("k"));
        vector<int> prior_type = as< vector<int> >(item.slot("prior.type"));
        vector<double> prior_1 = as< vector<double> >(item.slot("prior_1"));
        vector<double> prior_2 = as< vector<double> >(item.slot("prior_2"));
        NumericMatrix dat = item.slot("dat");
        NumericMatrix item_Q;
        if(itemclass == 10) item_Q = as<NumericMatrix>(item.slot("item.Q"));
        switch(itemclass){
            case -999: //custom group
                break;
            case -1 :
                _dgroupLCA(tmpgrad, tmphess, item, theta, estHess);
                break;
            case 0 :
                if(EMcomplete){
                    _dgroupEMCD(tmpgrad, tmphess, item, theta, estHess);
                } else {
                    _dgroupEM(tmpgrad, tmphess, item, theta, itemtrace, prior, estHess);
                }
                break;
            case 1 :
                d_dich(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, estHess);
                break;
            case 2 :
                d_poly(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat - 1, estHess);
                break;
            case 3 :
                has_mat = as<int>(item.slot("mat"));
                if(has_mat){
                    d_nominal2(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat, 0, estHess);
                } else {
                    d_nominal(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat, 0, estHess);
                }
                break;
            case 4 :
                d_nominal(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat, 0, estHess);
                break;
            case 6 :
                d_gpcmIRT(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat - 1, estHess);
                break;
            case 10 :
                d_lca(tmpgrad, tmphess, par, theta, item_Q, offterm(_,i), dat, N, nfact2, estHess);
                break;
            case 11 :
                d_ggum(tmpgrad, tmphess, par, theta, dat, N, nfact2, ncat, estHess);
                break;
            default :
            	d_numerical(tmpgrad, tmphess, par, theta, offterm(_,i), dat, N, nfact2, ncat, k, estHess, itemclass);
                break;
        }
        vector<int> parnum = as< vector<int> >(item.slot("parnum"));
        int where = parnum[0] - 1;
        for(int len = 0; len < par_size; ++len){
            if(useIprior)
                if(prior_type[len])
                    d_priors(tmpgrad, tmphess, len, prior_type[len], prior_1[len], prior_2[len], par[len]);
            grad[where + len] = tmpgrad[len];
            if(estHess){
                for(int len2 = 0; len2 < par_size; ++len2)
                    hess(where + len, where + len2) = tmphess(len, len2);
            }
        }
    }
}

RcppExport SEXP computeDPars(SEXP Rpars, SEXP RTheta, SEXP Roffterm,
    SEXP Rnpars, SEXP RestHess, SEXP RUSEFIXED, SEXP REM, SEXP REMcomplete)
{
    BEGIN_RCPP

    const List gpars(Rpars);
    const List gTheta(RTheta);
    const NumericMatrix offterm(Roffterm);
    const NumericMatrix dummy(1,1);
    const int nitems = offterm.ncol();
    const int npars = as<int>(Rnpars);
    const int EMcomplete = as<bool>(REMcomplete);
    const int estHess = as<int>(RestHess);
    const int USEFIXED = as<int>(RUSEFIXED);
    const int EM = as<int>(REM);
    vector<double> grad(npars);
    vector<double> dummy2(npars);
    int dim = 0; 
    if(estHess) dim = npars;
    NumericMatrix hess(dim, dim);

    for(int group = 0; group < gpars.length(); ++group){
        List pars = gpars[group];
        NumericMatrix Theta = gTheta[group];
        _computeDpars(grad, hess, pars, Theta, offterm, dummy, dummy2, nitems, npars,
            estHess, USEFIXED, EM, EMcomplete, true);
    }

    List ret;
    ret["grad"] = wrap(grad);
    ret["hess"] = hess;
    return(ret);

    END_RCPP
}

RcppExport SEXP computeInfo(SEXP Rpars, SEXP RTheta, SEXP RgPrior, SEXP Rgprior,
    SEXP RgPriorbetween, SEXP Rtabdata, SEXP Rrs, SEXP Rsitems, SEXP Ritemloc,
    SEXP Rgitemtrace, SEXP Rnpars, SEXP Rwmiss, SEXP Risbifactor, SEXP Riscross)
{
    BEGIN_RCPP

    List gpars(Rpars);
    const List gitemtrace(Rgitemtrace);
    const List gprior(Rgprior);
    const NumericMatrix gPrior(RgPrior); //cols are groups
    const NumericMatrix Theta(RTheta);
    const IntegerMatrix tabdata(Rtabdata);
    const IntegerMatrix sitems(Rsitems);
    const vector<int> itemloc = as< vector<int> >(Ritemloc);
    const NumericMatrix rs(Rrs); //group stacked
    const vector<double> wmiss = as< vector<double> >(Rwmiss);
    const NumericMatrix gPriorbetween(RgPriorbetween);
    const vector<double> vone(1.0, 1);
    const int nfact = Theta.ncol();
    const int nquad = Theta.nrow();
    const int J = itemloc[itemloc.size()-1] - 1;
    const int nitems = itemloc.size() - 1;
    const int npars = as<int>(Rnpars);
    const int isbifactor = as<int>(Risbifactor);
    const int iscross = as<int>(Riscross);
    const int ngroups = gpars.length();
    const int npat = tabdata.nrow();
    const bool Etable = true;
    NumericMatrix prior = gprior[0];
    const int nsfact = prior.ncol();
    const int nsquad = prior.nrow();
    const int npquad = gPriorbetween.ncol();
    IntegerMatrix dat(1, J);
    NumericMatrix Igrad(npars, npars), IgradP(npars, npars), offterm(1, nitems);

    for(int pat = 0; pat < npat; ++pat){
        for(int i = 0; i < J; ++i)
            dat(0, i) = tabdata(pat, i);
        for(int g = 0; g < ngroups; ++g){
            NumericMatrix itemtrace = gitemtrace[g];
            NumericVector tmpvec = gPrior(_,g);
            vector<double> Prior = as< vector<double> >(tmpvec);
            vector<double> expected(1), r1g(nquad), 
                r1vec(nquad*J), r2vec(npquad), r3vec(nsquad*nsfact);
            if(isbifactor){
                NumericMatrix prior = gprior[g];
                NumericVector tmpvec = gPriorbetween(g,_);
                vector<double> Priorbetween = as< vector<double> >(tmpvec);
               _Estepbfactor(expected, r1vec, r2vec, r3vec, itemtrace, prior, Priorbetween, vone,
                    dat, sitems, wmiss, Etable);
            } else {
                _Estep(expected, r1vec, r1g, Prior, vone, dat, itemtrace, wmiss, Etable);
            }
            NumericMatrix r1 = vec2mat(r1vec, nquad, J);
            List pars = gpars[g];
            if(iscross){
                for(int i = 0; i < nitems; ++i){
                    S4 item = pars[i];
                    NumericMatrix tmpmat(nquad, itemloc[i+1] - itemloc[i]);
                    for(int j = 0; j < tmpmat.ncol(); ++j)
                        for(int n = 0; n < nquad; ++n)
                            tmpmat(n,j) = r1(n, itemloc[i] + j - 1);
                    item.slot("dat") = tmpmat;
                    pars[i] = item;
                }
                S4 item = pars[nitems];
                item.slot("dat") = dat;
                item.slot("rr") = r1g;
                if(isbifactor){
                    NumericVector r2 = wrap(r2vec);
                    NumericMatrix r3 = vec2mat(r3vec, nsquad, nsfact);
                    item.slot("rrb") = r2;
                    item.slot("rrs") = r3;
                }
                pars[nitems] = item;
                NumericMatrix hess(npars, npars);
                vector<double> grad(npars);
                _computeDpars(grad, hess, pars, Theta, offterm, itemtrace, Prior,
                              nitems, npars, 0, 0, 1, true, false);
                add2outer(Igrad, grad, rs(g, pat));
            } else {
                for(int i = 0; i < nitems; ++i){
                    S4 item = pars[i];
                    NumericMatrix tmpmat(1, itemloc[i+1] - itemloc[i]);
                    for(int j = 0; j < tmpmat.ncol(); ++j){
                        tmpmat(0,j) = dat(0, itemloc[i] + j - 1);
                    }
                    item.slot("dat") = tmpmat;
                    pars[i] = item;
                }
                S4 item = pars[nitems];
                item.slot("dat") = dat;
                item.slot("rr") = wrap(1.0);
                // if(isbifactor){


                // }
                pars[nitems] = item;
                vector<double> w(nquad), grad(npars);
                for(int i = 0; i < J; ++i){
                    if(dat(0, i)){
                        for(int n = 0; n < nquad; ++n)
                            w[n] = r1(n,i);
                        break;
                    }
                }
                NumericMatrix theta(1, nfact);
                for(int n = 0; n < nquad; ++n){
                    for(int i = 0; i < nfact; ++i)
                        theta(0,i) = Theta(n,i);
                    NumericMatrix hess(npars, npars);
                    vector<double> tmpgrad(npars);
                    _computeDpars(tmpgrad, hess, pars, theta, offterm, itemtrace, Prior,
                                  nitems, npars, 0, 0, 1, true, false);
                    add2outer(IgradP, tmpgrad, rs(g,pat) * w[n]);
                    for(int j = 0; j < npars; ++j)
                        grad[j] += tmpgrad[j] * w[n];
                }
                add2outer(Igrad, grad, rs(g,pat));
            }
        }
    }

    List ret;
    ret["Igrad"] = Igrad;
    ret["IgradP"] = IgradP;
    return(ret);
    END_RCPP
}

RcppExport SEXP computeGradient(SEXP Rpars, SEXP RTheta, SEXP RgPrior,
    SEXP Rgprior, SEXP RgPriorbetween, SEXP Rtabdata, SEXP Rsitems,
    SEXP Ritemloc, SEXP Rgitemtrace, SEXP Rwmiss, SEXP Rnpars, SEXP Risbifactor)
{
    BEGIN_RCPP

    List gpars(Rpars);
    const List gitemtrace(Rgitemtrace);
    const List gprior(Rgprior);
    const NumericMatrix gPrior(RgPrior); //cols are groups
    const NumericMatrix Theta(RTheta);
    const IntegerMatrix tabdata(Rtabdata);
    const IntegerMatrix sitems(Rsitems);
    const vector<int> itemloc = as< vector<int> >(Ritemloc);
    const NumericMatrix gPriorbetween(RgPriorbetween);
    const vector<double> wmiss = as< vector<double> >(Rwmiss);
    const vector<double> vone(1.0, 1);
    const int nquad = Theta.nrow();
    const int J = itemloc[itemloc.size()-1] - 1;
    const int nitems = itemloc.size() - 1;
    const int npars = as<int>(Rnpars);
    const int isbifactor = as<int>(Risbifactor);
    const int ngroups = gpars.length();
    const int npat = tabdata.nrow();
    const bool Etable = true;
    NumericMatrix prior = gprior[0];
    const int nsfact = prior.ncol();
    const int nsquad = prior.nrow();
    const int npquad = gPriorbetween.ncol();
    IntegerMatrix dat(1, J);
    NumericMatrix offterm(1, nitems), gradient(npat, npars);

    for(int pat = 0; pat < npat; ++pat){
        for(int i = 0; i < J; ++i)
            dat(0, i) = tabdata(pat, i);
        for(int g = 0; g < ngroups; ++g){
            NumericMatrix itemtrace = gitemtrace[g];
            NumericVector tmpvec = gPrior(_,g);
            vector<double> Prior = as< vector<double> >(tmpvec);
            vector<double> expected(1), r1g(nquad),
                r1vec(nquad*J), r2vec(npquad), r3vec(nsquad*nsfact);
            if(isbifactor){
                NumericMatrix prior = gprior[g];
                NumericVector tmpvec = gPriorbetween(g,_);
                vector<double> Priorbetween = as< vector<double> >(tmpvec);
               _Estepbfactor(expected, r1vec, r2vec, r3vec, itemtrace, prior, Priorbetween, vone,
                    dat, sitems, wmiss, Etable);
            } else {
                _Estep(expected, r1vec, r1g, Prior, vone, dat, itemtrace, wmiss, Etable);
            }
            NumericMatrix r1 = vec2mat(r1vec, nquad, J);
            List pars = gpars[g];
            for(int i = 0; i < nitems; ++i){
                S4 item = pars[i];
                NumericMatrix tmpmat(nquad, itemloc[i+1] - itemloc[i]);
                for(int j = 0; j < tmpmat.ncol(); ++j)
                    for(int n = 0; n < nquad; ++n)
                        tmpmat(n,j) = r1(n, itemloc[i] + j - 1);
                item.slot("dat") = tmpmat;
                pars[i] = item;
            }
            S4 item = pars[nitems];
            item.slot("dat") = dat;
            item.slot("rr") = r1g;
            if(isbifactor){
                NumericVector r2 = wrap(r2vec);
                NumericMatrix r3 = vec2mat(r3vec, nsquad, nsfact);
                item.slot("rrb") = r2;
                item.slot("rrs") = r3;
            }
            pars[nitems] = item;
            NumericMatrix hess(npars, npars);
            vector<double> grad(npars);
            _computeDpars(grad, hess, pars, Theta, offterm, itemtrace, Prior,
                          nitems, npars, 0, 0, 1, true,true);
            for(int i = 0; i < npars; ++i){
                gradient(pat, i) += grad[i];
            }
        }
    }

    return(gradient);
    END_RCPP
}
