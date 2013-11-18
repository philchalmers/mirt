#include"Misc.h"

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

void itemTrace(vector<double> &P, vector<double> &Pstar, const vector<double> &a, const double *d, 
        const NumericMatrix &Theta, const double *g, const double *u, const NumericVector &ot)
{	
    const int nquad = Theta.nrow();
    const int nfact = Theta.ncol();
    const int USEOT = ot.length() > 1;

	for (int i = 0; i <	nquad; ++i){
        double z = *d;        
    	for (int j = 0; j <	nfact; ++j)		
			z += a[j] * Theta(i,j); 
        if(USEOT) z += ot(i);
        if(z > ABS_MAX_Z) z = ABS_MAX_Z;
        else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
        Pstar[i] = 1.0 / (1.0 + exp(-z));
    	P[i] = *g + (*u - *g) * Pstar[i];
	}
}

RcppExport SEXP reloadPars(SEXP Rlongpars, SEXP Rpars, SEXP Rngroups, SEXP RJ)
{    
    BEGIN_RCPP
	const NumericVector longpars(Rlongpars);
    List pars(Rpars);
    const int ngroups = as<int>(Rngroups);
    const int J = as<int>(RJ);
    int ind = 0;

    for(int g = 0; g < ngroups; ++g){
        List glist = pars[g];
        for(int i = 0; i < (J+1); ++i){
            S4 item = glist[i];
            NumericVector p = item.slot("par");
            int len = p.length();
            for(int j = 0; j < len; ++j)
                p(j) = longpars(ind+j);
            ind += len;
            item.slot("par") = p;
            glist[i] = item;
        }        
        pars[g] = glist;
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
    const NumericVector log_den0(Rlog_den0);
    const NumericVector log_den1(Rlog_den1);
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
        Sum0[i] = rs0 + log_den0(i);
        Sum1[i] = rs1 + log_den1(i);
    }
	
    ret["total_0"] = wrap(Sum0);
    ret["total_1"] = wrap(Sum1);
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

SEXP vec2mat(vector<double> &x, const int *nrow, const int *ncol) {
  NumericVector output = wrap(x);
  output.attr("dim") = Dimension(*nrow, *ncol);
  return(output);
}
