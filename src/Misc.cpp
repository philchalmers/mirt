#include "Misc.h"

// hack fix
void R_init_mirt(DllInfo* info) {
    R_registerRoutines(info, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}

NumericMatrix polyOuter(const NumericMatrix &Thetas, const vector<double> &Pk,
	const vector<double> &Pk_1, const vector<double> &PQ_1, const vector<double> &PQ,
	const vector<double> &dif1sq, const vector<double> &dif1)
{
	const int nfact = Thetas.ncol();
	NumericMatrix d2Louter(nfact,nfact), outer(nfact,nfact);
	vector<double> temp(nfact);

	for(int n = 0; n < Thetas.nrow(); ++n){
		for(int i = 0; i < nfact; ++i)
			for(int j = 0; j < nfact; ++j)
				outer(i,j) = Thetas(n,i) * Thetas(n,j);
		for(int i = 0; i < nfact; ++i)
			temp[i] =  (PQ_1[n] * Thetas(n,i) - PQ[n] * Thetas(n,i));
		for(int i = 0; i < nfact; ++i)
			for(int j = 0; j < nfact; ++j)
				d2Louter(i,j) += (-1) * dif1sq[n] * temp[i] * temp[j] +
				    (dif1[n] * (Pk_1[n] * (1.0 - Pk_1[n]) * (1.0 - 2.0 * Pk_1[n]) *
				    outer(i,j) - Pk[n] * (1.0 - Pk[n]) * (1.0 - 2.0 * Pk[n]) * outer(i,j)));
	}
	return d2Louter;
}

RcppExport SEXP reloadPars(SEXP Rlongpars, SEXP Rpars, SEXP Rngroups, SEXP RJ,
                           SEXP Rnclasspars)
{
    BEGIN_RCPP
	const NumericVector longpars(Rlongpars);
    List pars(Rpars);
    const int ngroups = as<int>(Rngroups);
    const int J = as<int>(RJ);
    const vector<int> nclasspars = as< vector<int> >(Rnclasspars);
    int ind = 0;

    for(int g = 0; g < ngroups; ++g){
        List glist = pars[g];
        for(int i = 0; i < (J+1); ++i){
            S4 item = glist[i];
            NumericVector par(nclasspars[i]);
            for(int j = 0; j < nclasspars[i]; ++j){
                par(j) = longpars(ind);
                ++ind;
            }
            item.slot("par") = par;
        }
    }

    return(pars);
	END_RCPP
}

RcppExport SEXP denRowSums(SEXP Rfulldata, SEXP Ritemtrace0, SEXP Ritemtrace1,
    SEXP Rlog_den0, SEXP Rlog_den1)
{
    BEGIN_RCPP

    const IntegerMatrix fulldata(Rfulldata);
    const NumericMatrix itemtrace0(Ritemtrace0);
    const NumericMatrix itemtrace1(Ritemtrace1);
    const vector<double> log_den0 = as< vector<double> >(Rlog_den0);
    const vector<double> log_den1 = as< vector<double> >(Rlog_den1);
    List ret;
    vector<double> Sum0(fulldata.nrow()), Sum1(fulldata.nrow());

    for(int i = 0; i < fulldata.nrow(); ++i){
        double rs0 = 0.0;
        double rs1 = 0.0;
        for(int j = 0; j < fulldata.ncol(); ++j){
            if(fulldata(i,j)){
                rs0 += log(itemtrace0(i,j));
                rs1 += log(itemtrace1(i,j));
            }
        }
        Sum0[i] = rs0 + log_den0[i];
        Sum1[i] = rs1 + log_den1[i];
    }

    ret["total_0"] = wrap(Sum0);
    ret["total_1"] = wrap(Sum1);
    return(ret);
	END_RCPP
}

RcppExport SEXP sumExpected(SEXP Rtdata, SEXP Rtabdata, SEXP Rrwmeans, SEXP Rnitems)
{
    BEGIN_RCPP

    const IntegerMatrix tdata(Rtdata);
    const IntegerMatrix tabdata(Rtabdata);
    const NumericVector rwmeans(Rrwmeans);
    const int nitems = as<int>(Rnitems);
    const int N = tdata.ncol();
    const int n = tabdata.nrow();
    const int J = tdata.nrow();
    vector<double> expected(n);

    for(int i = 0; i < n; ++i){
        int count = 0;
        double tempexp = 0.0;
        for(int NN = 0; NN < N; ++NN){
            int tmp = 0;
            for(int j = 0; j < J; ++j)
                tmp += tabdata(i, j) == tdata(j, NN);
            if(tmp == nitems){
                count += 1;
                tempexp = tempexp + rwmeans(NN);
            }
        }
        if(count) expected[i] = tempexp / count;
    }

    return(wrap(expected));

    END_RCPP
}

RcppExport SEXP buildXi2els(SEXP Rdim1, SEXP Rdim2, SEXP Rnitems, SEXP REIs,
    SEXP REIs2, SEXP RPrior)
{
    BEGIN_RCPP
    const int dim1 = as<int>(Rdim1);
    const int dim2 = as<int>(Rdim2);
    const int nitems = as<int>(Rnitems);
    const NumericMatrix EIs(REIs);
    const NumericMatrix EIs2(REIs2);
    const vector<double> Prior = as< vector<double> >(RPrior);
    const int N = EIs.nrow();
    NumericMatrix Xi11(dim1, dim1);
    NumericMatrix Xi12(dim1, dim2);
    NumericMatrix Xi22(dim2, dim2);

    for(int i = 0; i < nitems; ++i){
        for(int j = 0; j < nitems; ++j){
            if(i >= j){
                double pa = 0, pb = 0, pab = 0;
                for(int n = 0; n < N; ++n){
                    pa += EIs(n,i) * Prior[n];
                    pb += EIs(n,j) * Prior[n];
                }
                if(i == j){
                    for(int n = 0; n < N; ++n)
                        pab += EIs2(n,i) * Prior[n];
                } else {
                    for(int n = 0; n < N; ++n)
                        pab += EIs(n,i) * EIs(n,j) * Prior[n];
                }
                Xi11(i,j) = pab - pa*pb;
                Xi11(j,i) = Xi11(i,j);
            }
        }
    }
    for(int k = 0; k < nitems; ++k){
        int ind = 0;
        for(int i = 0; i < nitems; ++i){
            for(int j = 0; j < nitems; ++j){
                if(i < j){
                    double pab = 0, pc = 0, pabc = 0;
                    for(int n = 0; n < N; ++n){
                        pab += EIs(n,i) * EIs(n,j) * Prior[n];
                        pc += EIs(n,k) * Prior[n];
                    }
                    if(i == k){
                        for(int n = 0; n < N; ++n)
                            pabc += EIs2(n,i) * EIs(n,j) * Prior[n];
                    } else if(j == k){
                        for(int n = 0; n < N; ++n)
                            pabc += EIs(n,i) * EIs2(n,j) * Prior[n];
                    } else {
                        for(int n = 0; n < N; ++n)
                            pabc += EIs(n,i) * EIs(n,j) * EIs(n,k) * Prior[n];
                    }
                    Xi12(k,ind) = pabc - pab*pc;
                    ++ind;
                }
            }
        }
    }
    int ind1 = 0;
    for(int k = 0; k < nitems; ++k){
        for(int l = 0; l < nitems; ++l){
            if(k < l){
                int ind2 = 0;
                for(int i = 0; i < nitems; ++i){
                    for(int j = 0; j < nitems; ++j){
                        if(i < j){
                            double pab = 0, pcd = 0, pabcd = 0;
                            for(int n = 0; n < N; ++n){
                                pab += EIs(n,i) * EIs(n,j) * Prior[n];
                                pcd += EIs(n,k) * EIs(n,l) * Prior[n];
                            }
                            if((i == k) & (j == l)){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs2(n,j) * Prior[n];
                            } else if(i == k){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs(n,j) * EIs(n,l) * Prior[n];
                            } else if(j == k){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs2(n,j) * EIs(n,l) * Prior[n];
                            } else if(i == l){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs(n,j) * EIs(n,k) * Prior[n];
                            } else if(j == l){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs2(n,j) * EIs(n,k) * Prior[n];
                            } else {
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs(n,j) * EIs(n,k) * EIs(n,l) * Prior[n];
                            }
                            Xi22(ind1, ind2) = pabcd - pab*pcd;
                            ++ind2;
                        }
                    }
                }
                ++ind1;
            }
        }
    }

    List ret;
    ret["Xi11"] = Xi11;
    ret["Xi12"] = Xi12;
    ret["Xi22"] = Xi22;
    return(ret);

    END_RCPP
}

RcppExport SEXP buildXi2els_C2(SEXP Rdim1, SEXP Rdim2, SEXP Rnitems0, 
	SEXP Rnitems, SEXP RPIs, SEXP REIs, SEXP REIs2, SEXP RPrior, 
	SEXP Rabcats, SEXP Rabcats2)
{
    BEGIN_RCPP
    const int dim1 = as<int>(Rdim1);
    const int dim2 = as<int>(Rdim2);
    const int nitems0 = as<int>(Rnitems0);
    const int nitems = as<int>(Rnitems);
    const NumericMatrix PIs(RPIs);
    const NumericMatrix EIs(REIs);
    const NumericMatrix EIs2(REIs2);
    const vector<double> Prior = as< vector<double> >(RPrior);
    const vector<int> abcats = as< vector<int> >(Rabcats);
    const vector<int> abcats2 = as< vector<int> >(Rabcats2);
    const int N = EIs.nrow();
    NumericMatrix Xi11(dim1, dim1);
    NumericMatrix Xi12(dim1, dim2);
    NumericMatrix Xi22(dim2, dim2);

    for(int i = 0; i < nitems0; ++i){
        for(int j = 0; j < nitems0; ++j){
            if(i >= j){
                double pa = 0, pb = 0, pab = 0;
                for(int n = 0; n < N; ++n){
                    pa += PIs(n,i) * Prior[n];
                    pb += PIs(n,j) * Prior[n];
                }
                if(abcats[i] == abcats[j]){
                	if(abcats2[i] == abcats2[j])
	                    for(int n = 0; n < N; ++n)
	                        pab += PIs(n,i) * Prior[n];
                } else {
                    for(int n = 0; n < N; ++n)
                		pab += PIs(n,i) * PIs(n,j) * Prior[n];
                }
                Xi11(i,j) = pab - pa*pb;
                Xi11(j,i) = Xi11(i,j);
            }
        }
    }
    for(int k = 0; k < nitems0; ++k){
    	int ind = 0;
        for(int i = 0; i < nitems; ++i){
            for(int j = 0; j < nitems; ++j){
                if(i < j){
                    double pab = 0, pc = 0, pabc = 0;
                    for(int n = 0; n < N; ++n){
                        pab += EIs(n,i) * EIs(n,j) * Prior[n];
                        pc += PIs(n,k) * Prior[n];
                    }
                    if(i == abcats[k]){ 
	                    for(int n = 0; n < N; ++n)
	                    	pabc += abcats2[k] * PIs(n,k) * EIs(n,j) * Prior[n];
                    } else if(j == abcats[k]){
	                    for(int n = 0; n < N; ++n)
	                        pabc += abcats2[k] * PIs(n,k) * EIs(n,i) * Prior[n];
                    } else {
                        for(int n = 0; n < N; ++n)
                    		pabc += EIs(n,i) * EIs(n,j) * PIs(n,k) * Prior[n];
                    }
                    Xi12(k,ind) = pabc - pab*pc;
                    ++ind;
                }
            }
        }
    }
    int ind1 = 0;
    for(int k = 0; k < nitems; ++k){
        for(int l = 0; l < nitems; ++l){
            if(k < l){
                int ind2 = 0;
                for(int i = 0; i < nitems; ++i){
                    for(int j = 0; j < nitems; ++j){
                        if(i < j){
                            double pab = 0, pcd = 0, pabcd = 0;
                            for(int n = 0; n < N; ++n){
                                pab += EIs(n,i) * EIs(n,j) * Prior[n];
                                pcd += EIs(n,k) * EIs(n,l) * Prior[n];
                            }
                            if((i == k) & (j == l)){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs2(n,j) * Prior[n];
                            } else if(i == k){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs(n,j) * EIs(n,l) * Prior[n];
                            } else if(j == k){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs2(n,j) * EIs(n,l) * Prior[n];
                            } else if(i == l){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs2(n,i) * EIs(n,j) * EIs(n,k) * Prior[n];
                            } else if(j == l){
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs2(n,j) * EIs(n,k) * Prior[n];
                            } else {
                                for(int n = 0; n < N; ++n)
                                    pabcd += EIs(n,i) * EIs(n,j) * EIs(n,k) * EIs(n,l) * Prior[n];
                            }
                            Xi22(ind1, ind2) = pabcd - pab*pcd;
                            ++ind2;
                        }
                    }
                }
                ++ind1;
            }
        }
    }

    List ret;
    ret["Xi11"] = Xi11;
    ret["Xi12"] = Xi12;
    ret["Xi22"] = Xi22;
    return(ret);

    END_RCPP
}

double antilogit(const double *x){
    double ret;
    if(*x > 998.0) ret = 1.0;
    else if(*x < -998.0) ret = 0.0;
    else ret = 1.0 / (1.0 + exp(-1.0 * (*x)));
    return(ret);
}

double vecsum(const vector<double> &x)
{
    double sum = 0.0;
    const int size = x.size();
    for(int i = 0; i < size; ++i)
        sum += x[i];
    return(sum);
}

SEXP vec2mat(vector<double> &x, const int &nrow, const int &ncol) {
  NumericVector output = wrap(x);
  output.attr("dim") = Dimension(nrow, ncol);
  return(output);
}

RcppExport SEXP respSample(SEXP Ritemtrace)
{
    BEGIN_RCPP

    const NumericMatrix itemtrace(Ritemtrace);
    const int ncat = itemtrace.ncol();
    const int N = itemtrace.nrow();
    const NumericVector unif = runif(N);
    std::vector<int> resp(N);

    for(int i = 0; i < N; ++i){
        double csum = itemtrace(i, 0);
        int cat = 0;
        while(csum < unif(i)){
            ++cat;
            if(cat == ncat) break;
            csum += itemtrace(i, cat);
        }
        resp[i] = cat;
    }

    return(wrap(resp));
    END_RCPP
}
