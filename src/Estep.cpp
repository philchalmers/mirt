#include"Misc.h"

//Estep for mirt
RcppExport SEXP Estep(SEXP Ritemtrace, SEXP Rprior, SEXP RX, SEXP Rr)
{
    BEGIN_RCPP

    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<int> r = as< vector<int> >(Rr);
    const IntegerMatrix data(RX);
    const NumericMatrix itemtrace(Ritemtrace);
    const int nquad = prior.size();
    const int nitems = data.ncol();
    const int npat = r.size();
    vector<double> expected(npat, 0.0);
    vector<double> r1vec(nquad*nitems, 0.0);
    List ret;

    // Begin main function body
    for (int pat = 0; pat < npat; ++pat){
        vector<double> posterior(nquad,1.0);
        for(int q = 0; q < nquad; ++q)
            posterior[q] = posterior[q] * prior[q];
        for (int item = 0; item < nitems; ++item)
            if(data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    posterior[q] = posterior[q] * itemtrace(q,item);
        double expd = 0.0;
        for(int i = 0; i < nquad; ++i)
            expd += posterior[i];
        expected[pat] = expd;
        for(int q = 0; q < nquad; ++q)
            posterior[q] = r[pat] * posterior[q] / expd;
        for (int item = 0; item < nitems; ++item)
            if (data(pat,item))
                for(int q = 0; q < nquad; ++q)
                    r1vec[q + item*nquad] += posterior[q];
    } //end main

    NumericMatrix r1 = vec2mat(r1vec, nquad, nitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    return(ret);

    END_RCPP
}


//Estep for bfactor
RcppExport SEXP Estepbfactor(SEXP Ritemtrace, SEXP Rprior, SEXP RPriorbetween, SEXP RX,
    SEXP Rr, SEXP Rsitems)
{
    BEGIN_RCPP

    List ret;
    const NumericMatrix itemtrace(Ritemtrace);
    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> Priorbetween = as< vector<double> >(RPriorbetween);
    const vector<int> r = as< vector<int> >(Rr);
    const IntegerMatrix data(RX);
    const IntegerMatrix sitems(Rsitems);
    const int sfact = sitems.ncol();
    const int nitems = data.ncol();
    const int npquad = prior.size();
    const int nbquad = Priorbetween.size();
    const int nquad = nbquad * npquad;
    const int npat = r.size();

    vector<double> expected(npat);
    vector<double> r1vec(nquad*nitems*sfact, 0.0);

    // Begin main function body here
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
        for (int item = 0; item < nitems; ++item)
            if (data(pat,item))
                for (int fact = 0; fact < sfact; ++fact)
                    for(int q = 0; q < nquad; ++q)
                        r1vec[q + fact*nquad*nitems + nquad*item] += posterior[q + fact*nquad];
    }   //end main

    int nsitems = sfact * nitems;
    NumericMatrix r1 = vec2mat(r1vec, nquad, nsitems);
    ret["r1"] = r1;
    ret["expected"] = wrap(expected);
    return(ret);

    END_RCPP
}

//EAP estimates used in multipleGroup
RcppExport SEXP EAPgroup(SEXP Ritemtrace, SEXP Rtabdata, SEXP RTheta, SEXP Rprior, SEXP Rmu)
{
    BEGIN_RCPP

    const NumericMatrix itemtrace(Ritemtrace);
    const NumericMatrix tabdata(Rtabdata);
    const NumericMatrix Theta(RTheta);
    const vector<double> prior = as< vector<double> >(Rprior);
    const vector<double> mu = as< vector<double> >(Rmu);
    const int n = prior.size(); //nquad
    const int nitems = tabdata.ncol();
    const int nfact = Theta.ncol();
    NumericMatrix scores(tabdata.nrow(), nfact), scores2(tabdata.nrow(), nfact*(nfact + 1)/2);

    for(int pat = 0; pat < tabdata.nrow(); ++pat){

        vector<double> L(n, 1.0);
        for(int j = 0; j < n; ++j){
            for(int i = 0; i < nitems; ++i)
                if(tabdata(pat, i))
                    L[j] = L[j] * itemtrace(j, i);
        }

        vector<double> thetas(nfact, 0.0);
        double denom = 0.0;
        for(int j = 0; j < n; ++j)
            denom += (L[j] * prior[j]);

        for(int k = 0; k < nfact; ++k){
            for(int j = 0; j < n; ++j)
                thetas[k] += Theta(j, k) * L[j] * prior[j] / denom;
            scores(pat, k) = thetas[k];
        }

        int ind = 0;
        vector<double> thetas2(nfact*(nfact+1)/2, 0.0);
        for(int i = 0; i < nfact; ++i){
            for(int k = 0; k < nfact; ++k){
                if(i <= k){
                    for(int j = 0; j < n; ++j)
                        thetas2[ind] += (Theta(j, i) - thetas[i]) * (Theta(j, k) - thetas[k]) *
                            L[j] * prior[j] / denom;
                    thetas2[ind] += (thetas[i] - mu[i]) * (thetas[k] - mu[k]);
                    scores2(pat, ind) = thetas2[ind];
                    ind += 1;
                }
            }
        }
    }

    List ret;
    ret["scores"] = scores;
    ret["scores2"] = scores2;
    return(ret);

    END_RCPP
}
