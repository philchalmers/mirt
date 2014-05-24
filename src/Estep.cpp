#include "Misc.h"
#include "Estep.h"

void _Estep(vector<double> &expected, vector<double> &r1vec, const vector<double> &prior,
    const vector<int> &r, const IntegerMatrix &data, const NumericMatrix &itemtrace,
    const int &ncores)
{
    const int nquad = prior.size();
    const int nitems = data.ncol();
    const int npat = r.size();
    #ifdef SUPPORT_OPENMP
    if(nquad * nitems > 1000){
        omp_set_num_threads(ncores);
    } else {
        omp_set_num_threads(1);
    }
    #endif

//    #pragma omp parallel for
    for (int pat = 0; pat < npat; ++pat){
        vector<double> posterior(nquad,1.0);
        for(int q = 0; q < nquad; ++q)
            posterior[q] = posterior[q] * prior[q];
        for (int item = 0; item < nitems; ++item)
            if(data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    posterior[q] *= itemtrace(q,item);
        double expd = 0.0;
        for(int i = 0; i < nquad; ++i)
            expd += posterior[i];
        expected[pat] = expd;
        for(int q = 0; q < nquad; ++q)
            posterior[q] = r[pat] * posterior[q] / expd;
       #pragma omp critical
        for (int item = 0; item < nitems; ++item)
            if (data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    r1vec[q + item*nquad] += posterior[q];
    } //end main

}

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX, SEXP Rr, SEXP Rncores)
{
    BEGIN_RCPP

    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<int> r = as< vector<int> >(Rr);
    const IntegerMatrix data(RX);
    const NumericMatrix itemtrace(Ritemtrace);
    const int ncores = as<int>(Rncores);
    const int nquad = prior.size();
    const int nitems = data.ncol();
    const int npat = r.size();
    vector<double> expected(npat, 0.0);
    vector<double> r1vec(nquad*nitems, 0.0);
    List ret;

    _Estep(expected, r1vec, prior, r, data, itemtrace, ncores);
    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    return(ret);

    END_RCPP
}

void _Estepbfactor(vector<double> &expected, vector<double> &r1, vector<double> &ri,
    const NumericMatrix &itemtrace, const vector<double> &prior, const vector<double> &Priorbetween, 
    const vector<int> &r, const int &ncores, const IntegerMatrix &data, const IntegerMatrix &sitems,
    const vector<double> &Prior)
{
     #ifdef SUPPORT_OPENMP
    omp_set_num_threads(ncores);
    #endif
    const int sfact = sitems.ncol();
    const int nitems = data.ncol();
    const int npquad = prior.size();
    const int nbquad = Priorbetween.size();
    const int nquad = nbquad * npquad;
    const int npat = r.size();
    vector<double> r1vec(nquad*nitems*sfact, 0.0);

//#pragma omp parallel for
    for (int pat = 0; pat < npat; ++pat){
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
            vector<double> tempsum(nbquad, 0.0);
            for (int i = 0; i < npquad; ++i)
                for (int q = 0; q < nbquad; ++q)
                    tempsum[q] += L[q + i*nbquad];
            for (int i = 0; i < nbquad; ++i)
                Plk[i + fact*nbquad] = tempsum[i];
        }
        vector<double> Pls(nbquad, 1.0);
        for (int i = 0; i < nbquad; ++i){
            for(int fact = 0; fact < sfact; ++fact)
                Pls[i] = Pls[i] * Plk[i + fact*nbquad];
            expected[pat] += Pls[i] * Priorbetween[i];
        }
        for (int fact = 0; fact < sfact; ++fact)
            for (int i = 0; i < nbquad; ++i)
                Elk[i + fact*nbquad] = Pls[i] / Plk[i + fact*nbquad];
        for (int fact = 0; fact < sfact; ++fact)
            for (int i = 0; i < nquad; ++i)
                posterior[i + nquad*fact] = likelihoods[i + nquad*fact] * r[pat] * Elk[i % nbquad + fact*nbquad] /
                                            expected[pat];
        #pragma omp critical
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
    SEXP Rr, SEXP Rsitems, SEXP RPrior, SEXP Rncores)
{
    BEGIN_RCPP

    List ret;
    const NumericMatrix itemtrace(Ritemtrace);
    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> Priorbetween = as< vector<double> >(RPriorbetween);
    const vector<double> Prior = as< vector<double> >(RPrior);
    const vector<int> r = as< vector<int> >(Rr);
    const int ncores = as<int>(Rncores);
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

    _Estepbfactor(expected, r1vec, r2vec, itemtrace, prior, Priorbetween, r, ncores,
        data, sitems, Prior);
    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    ret["r2"] = wrap(r2vec);
    return(ret);

    END_RCPP
}

//EAP estimates used in multipleGroup
RcppExport SEXP EAPgroup(SEXP Rlogitemtrace, SEXP Rtabdata, SEXP RTheta, SEXP Rprior, SEXP Rmu,
    SEXP Rfull, SEXP Rr, SEXP Rncores)
{
    BEGIN_RCPP

    const int ncores = as<int>(Rncores);
    const int full = as<int>(Rfull);
    const vector<int> r = as< vector<int> >(Rr);
    #ifdef SUPPORT_OPENMP
    omp_set_num_threads(ncores);
    #endif
    const NumericMatrix logitemtrace(Rlogitemtrace);
    const IntegerMatrix tabdata(Rtabdata);
    const NumericMatrix Theta(RTheta);
    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> mu = as< vector<double> >(Rmu);
    const int n = prior.size(); //nquad
    const int N = tabdata.nrow();
    const int nitems = tabdata.ncol();
    const int nfact = Theta.ncol();
    vector<double> scores(N * nfact);
    vector<double> scores2(N * nfact*(nfact + 1)/2);

#pragma omp parallel for
    for(int pat = 0; pat < N; ++pat){

        vector<double> L(n, 0.0);
        for(int j = 0; j < n; ++j){
            for(int i = 0; i < nitems; ++i)
                if(tabdata(pat, i))
                    L[j] = L[j] + logitemtrace(j, i);
        }

        vector<double> thetas(nfact, 0.0);
        double denom = 0.0;
        for(int j = 0; j < n; ++j){
            L[j] = exp(L[j] + log(prior[j]));
            denom += L[j];
        }

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
                    thetas2[ind] += (thetas[i] - mu[i]) * (thetas[k] - mu[k]);
                    scores2[pat + ind*N] = thetas2[ind];
                    ind += 1;
                }
            }
        }
    }

    if(full){
        const int NN = std::accumulate(r.begin(), r.end(), 0);
        NumericMatrix fullscores(NN, nfact);
        NumericMatrix scoresmat = vec2mat(scores, N, nfact);
        int ind = 0;
        for(int pat = 0; pat < N; ++pat){
            for(int j = 0; j < r[pat]; ++j){
                for(int i = 0; i < nfact; ++i)
                    fullscores(ind, i) = scoresmat(pat, i);
            ++ind;
            }
        }
        return(fullscores);
    }
    List ret;
    ret["scores"] = vec2mat(scores, N, nfact);
    ret["scores2"] = vec2mat(scores2, N, nfact*(nfact + 1)/2);
    return(ret);

    END_RCPP
}
