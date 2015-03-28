#include "Misc.h"
#include "Estep.h"

const double ABSMIN = std::numeric_limits<double>::min();

void _Estep(vector<double> &expected, vector<double> &r1vec, const vector<double> &prior,
    const vector<double> &r, const IntegerMatrix &data, const NumericMatrix &itemtrace)
{
    const int nquad = prior.size();
    const int nitems = data.ncol();
    const int npat = r.size();

    for (int pat = 0; pat < npat; ++pat){
        if(r[pat] < 1e-10) continue;
        vector<double> posterior(nquad,1.0);
        for(int q = 0; q < nquad; ++q)
            posterior[q] = posterior[q] * prior[q];
        for (int item = 0; item < nitems; ++item)
            if(data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    posterior[q] *= itemtrace(q,item);
        const double maxp = *std::max_element(posterior.begin(), posterior.end());
        double expd = 0.0;
        for(int i = 0; i < nquad; ++i)
            expd += posterior[i]/maxp;
        expd *= maxp;
        if(expd > ABSMIN){
            for(int q = 0; q < nquad; ++q)
                posterior[q] = r[pat] * posterior[q] / expd;
        } else expd = ABSMIN;
        expected[pat] = expd;
        for (int item = 0; item < nitems; ++item)
            if (data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    r1vec[q + item*nquad] += posterior[q];
    } //end main

}

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX, SEXP Rr)
{
    BEGIN_RCPP

    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> r = as< vector<double> >(Rr);
    const IntegerMatrix data(RX);
    const NumericMatrix itemtrace(Ritemtrace);
    const int nquad = prior.size();
    const int nitems = data.ncol();
    const int npat = r.size();
    vector<double> expected(npat, 0.0);
    vector<double> r1vec(nquad*nitems, 0.0);
    List ret;

    _Estep(expected, r1vec, prior, r, data, itemtrace);
    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    return(ret);

    END_RCPP
}

void _Estep2(vector<double> &expected, vector<double> &r1vec, const NumericMatrix &prior,
    const IntegerMatrix &data, const NumericMatrix &itemtrace)
{
    const int nquad = prior.ncol();
    const int nitems = data.ncol();
    const int npat = data.nrow();

    for (int pat = 0; pat < npat; ++pat){
        vector<double> posterior(nquad,1.0);
        for(int q = 0; q < nquad; ++q)
            posterior[q] = posterior[q] * prior(pat,q);
        for (int item = 0; item < nitems; ++item)
            if(data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    posterior[q] *= itemtrace(q,item);
        const double maxp = *std::max_element(posterior.begin(), posterior.end());
        double expd = 0.0;
        for(int i = 0; i < nquad; ++i)
            expd += posterior[i]/maxp;
        expd *= maxp;
         if(expd > ABSMIN){
             for(int q = 0; q < nquad; ++q)
                 posterior[q] = posterior[q] / expd;
         } else expd = ABSMIN;
         expected[pat] = expd;
         for (int item = 0; item < nitems; ++item)
             if (data(pat,item))
                 for(int q = 0; q < nquad; ++q)
                     r1vec[q + item*nquad] += posterior[q];
     } //end main

}

//Estep for mirt
RcppExport SEXP Estep2(SEXP Ritemtrace, SEXP Rprior, SEXP RX)
{
    BEGIN_RCPP

    const NumericMatrix prior(Rprior);
    const IntegerMatrix data(RX);
    const NumericMatrix itemtrace(Ritemtrace);
    const int nquad = prior.ncol();
    const int nitems = data.ncol();
    const int npat = data.nrow();
    vector<double> expected(npat, 0.0);
    vector<double> r1vec(nquad*nitems, 0.0);
    List ret;

    _Estep2(expected, r1vec, prior, data, itemtrace);
    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    return(ret);

    END_RCPP
}

void _Estepbfactor(vector<double> &expected, vector<double> &r1, vector<double> &ri,
    const NumericMatrix &itemtrace, const vector<double> &prior, const vector<double> &Priorbetween,
    const vector<double> &r, const IntegerMatrix &data, const IntegerMatrix &sitems,
    const vector<double> &Prior)
{
    const int sfact = sitems.ncol();
    const int nitems = data.ncol();
    const int npquad = prior.size();
    const int nbquad = Priorbetween.size();
    const int nquad = nbquad * npquad;
    const int npat = r.size();
    vector<double> r1vec(nquad*nitems*sfact, 0.0);
    
    for (int pat = 0; pat < npat; ++pat){
        if(r[pat] < 1e-10) continue;
        vector<double> L(nquad), Elk(nbquad*sfact), posterior(nquad*sfact);
        vector<double> likelihoods(nquad*sfact, 1.0);
        for (int fact = 0; fact < sfact; ++fact){
            for (int item = 0; item < nitems; ++item){
                if (data(pat,item) && sitems(item,fact))
                    for (int k = 0; k < nquad; ++k)
                        likelihoods[k + nquad*fact] = likelihoods[k + nquad*fact] * itemtrace(k,item);
            }
        }
        vector<double> Plk(nbquad*sfact);
        for (int fact = 0; fact < sfact; ++fact){
            int k = 0;
            for (int q = 0; q < npquad; ++q){
                for (int i = 0; i < nbquad; ++i){
                    L[k] = likelihoods[k + nquad*fact] * prior[q];
                    ++k;
                }
            }
            const double maxL = *std::max_element(L.begin(), L.end());
            vector<double> tempsum(nbquad, 0.0);
            for (int i = 0; i < npquad; ++i)
                for (int q = 0; q < nbquad; ++q)
                    tempsum[q] += L[q + i*nbquad]/maxL;
            for (int i = 0; i < nbquad; ++i)
                Plk[i + fact*nbquad] = tempsum[i] * maxL;
        }
        vector<double> Pls(nbquad, 1.0);
        vector<double> PlsPlk(nbquad, 1.0);
        for (int i = 0; i < nbquad; ++i){
            for(int fact = 0; fact < sfact; ++fact)
                Pls[i] = Pls[i] * Plk[i + fact*nbquad];
            PlsPlk[i] = Pls[i] * Priorbetween[i];
        }
        double sumexp = 0.0;
        const double maxPlsPlk = *std::max_element(PlsPlk.begin(), PlsPlk.end());
        for (int i = 0; i < nbquad; ++i)
            sumexp += PlsPlk[i] / maxPlsPlk;
        expected[pat] = sumexp * maxPlsPlk;
        for (int fact = 0; fact < sfact; ++fact)
            for (int i = 0; i < nbquad; ++i)
                Elk[i + fact*nbquad] = Pls[i] / Plk[i + fact*nbquad];
        for (int fact = 0; fact < sfact; ++fact)
            for (int i = 0; i < nquad; ++i)
                posterior[i + nquad*fact] = likelihoods[i + nquad*fact] * r[pat] * Elk[i % nbquad + fact*nbquad] /
                                            expected[pat];
        for (int i = 0; i < nbquad; ++i)
            ri[i] += Pls[i] * r[pat] * Priorbetween[i] / expected[pat];
        for (int item = 0; item < nitems; ++item)
            if (data(pat,item))
                for (int fact = 0; fact < sfact; ++fact)
                    for(int q = 0; q < nquad; ++q)
                        r1vec[q + fact*nquad*nitems + nquad*item] += posterior[q + fact*nquad];
    }   //end main

    for (int item = 0; item < nitems; ++item)
        for (int fact = 0; fact < sfact; ++fact)
            if(sitems(item, fact))
                for(int q = 0; q < nquad; ++q)
                    r1[q + nquad*item] = r1vec[q + nquad*item + nquad*nitems*fact] * Prior[q];
}

//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP RPriorbetween, SEXP RX,
    SEXP Rr, SEXP Rsitems, SEXP RPrior)
{
    BEGIN_RCPP

    List ret;
    const NumericMatrix itemtrace(Ritemtrace);
    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> Priorbetween = as< vector<double> >(RPriorbetween);
    const vector<double> Prior = as< vector<double> >(RPrior);
    const vector<double> r = as< vector<double> >(Rr);
    const IntegerMatrix data(RX);
    const IntegerMatrix sitems(Rsitems);
    const int nitems = data.ncol();
    const int npquad = prior.size();
    const int nbquad = Priorbetween.size();
    const int nquad = nbquad * npquad;
    const int npat = r.size();
    vector<double> expected(npat);
    vector<double> r1vec(nquad*nitems, 0.0);
    vector<double> r2vec(nbquad, 0.0);

    _Estepbfactor(expected, r1vec, r2vec, itemtrace, prior, Priorbetween, r, 
        data, sitems, Prior);
    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    ret["r2"] = wrap(r2vec);
    return(ret);

    END_RCPP
}

RcppExport SEXP EAPgroup(SEXP Ritemtrace, SEXP Rtabdata, SEXP RTheta, SEXP Rprior, SEXP Rmu)
{
    BEGIN_RCPP
    
    const NumericMatrix itemtrace(Ritemtrace);
    const IntegerMatrix tabdata(Rtabdata);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix prior(Rprior);
    const NumericMatrix mu(Rmu);
    const int n = prior.ncol();
    const int N = tabdata.nrow();
    const int nitems = tabdata.ncol();
    const int nfact = Theta.ncol();
    vector<double> scores(N * nfact);
    vector<double> scores2(N * nfact*(nfact + 1)/2);

    for(int pat = 0; pat < N; ++pat){

        vector<double> L(n);
        for(int j = 0; j < n; ++j)
            L[j] = prior(pat, j);
        for(int j = 0; j < n; ++j){
            for(int i = 0; i < nitems; ++i)
                if(tabdata(pat, i))
                    L[j] *= itemtrace(j, i);
        }

        vector<double> thetas(nfact, 0.0);
        double denom = 0.0;
        const double maxL = *std::max_element(L.begin(), L.end());
        for(int j = 0; j < n; ++j)
            denom += L[j]/maxL;
        denom *= maxL;

        for(int k = 0; k < nfact; ++k){
            for(int j = 0; j < n; ++j)
                thetas[k] += Theta(j, k) * L[j] / denom;
            scores[pat + k*N] = thetas[k];
        }

        int ind = 0;
        vector<double> thetas2(nfact*(nfact+1)/2, 0.0);
        for(int i = 0; i < nfact; ++i){
            for(int k = 0; k < nfact; ++k){
                if(i <= k){
                    for(int j = 0; j < n; ++j)
                        thetas2[ind] += (Theta(j, i) - thetas[i]) * (Theta(j, k) - thetas[k]) *
                            L[j] / denom;
                    thetas2[ind] += (thetas[i] - mu(pat,i)) * (thetas[k] - mu(pat,k));
                    scores2[pat + ind*N] = thetas2[ind];
                    ind += 1;
                }
            }
        }
    }

    List ret;
    ret["scores"] = vec2mat(scores, N, nfact);
    ret["scores2"] = vec2mat(scores2, N, nfact*(nfact + 1)/2);
    return(ret);

    END_RCPP
}
