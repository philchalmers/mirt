#include"Misc.h"

const double ABS_MAX_Z = 30;

RcppExport SEXP dichOuter(SEXP RThetas, SEXP RPQ, SEXP RN)
{	
    BEGIN_RCPP
    NumericMatrix Thetas(RThetas);    
    NumericVector PQ(RPQ);
    NumericVector N(RN);
    const int nfact = Thetas.ncol();
	NumericMatrix ret(nfact,nfact);			

	for(int n = 0; n < N(0); ++n)
		for(int i = 0; i < nfact; ++i)
			for(int j = 0; j < nfact; ++j)
				ret(i,j) += Thetas(n,i) * Thetas(n,j) * PQ(n);
		
	return(ret);
	END_RCPP
}

NumericMatrix polyOuter(NumericMatrix Thetas, NumericVector Pk,
	NumericVector Pk_1, NumericVector PQ_1, NumericVector PQ, 
	NumericVector dif1sq, NumericVector dif1)
{
	const int nfact = Thetas.ncol();
	NumericMatrix d2Louter(nfact,nfact), outer(nfact,nfact);
	NumericVector temp(nfact);
	d2Louter.fill(0.0);
	
	for(int n = 0; n < Thetas.nrow(); ++n){
		for(int i = 0; i < nfact; ++i)
			for(int j = 0; j < nfact; ++j)
				outer(i,j) = Thetas(n,i) * Thetas(n,j);
		for(int i = 0; i < nfact; ++i)
			temp(i) =  (PQ_1(n) * Thetas(n,i) - PQ(n) * Thetas(n,i));
		for(int i = 0; i < nfact; ++i)
			for(int j = 0; j < nfact; ++j)				
				d2Louter(i,j) += (-1) * dif1sq(n) * temp(i) * temp(j) +  
				    (dif1(n) * (Pk_1(n) * (1.0 - Pk_1(n)) * (1.0 - 2.0 * Pk_1(n)) * 
				    outer(i,j) - Pk(n) * (1.0 - Pk(n)) * (1.0 - 2.0 * Pk(n)) * outer(i,j)));
	}
	return d2Louter;		
}

NumericVector itemTrace(NumericVector a, const double *d, 
        NumericMatrix Theta, const double *g, const double *u, NumericVector ot)
{	
    const int nquad = Theta.nrow();
    const int USEOT = ot.length() > 1;
	NumericVector P(nquad), z(nquad);
    z.fill(*d);

	for (int i = 0; i <	nquad; ++i){
		for (int j = 0; j <	Theta.ncol(); ++j)
			z(i) += a(j) * Theta(i,j);  
	}	
    if(USEOT){
        for (int i = 0; i < nquad; ++i)
            z(i) += ot(i);
    }
	for (int i = 0; i < nquad; ++i){
        if(z(i) > ABS_MAX_Z) z(i) = ABS_MAX_Z;
        else if(z(i) < -ABS_MAX_Z) z(i) = -ABS_MAX_Z;
		P(i) = *g + (*u - *g) / (1.0 + exp(-z(i)));
    }
	
	return P;		
}

RcppExport SEXP reloadPars(SEXP Rlongpars, SEXP Rpars, SEXP Rngroups, SEXP RJ)
{    
    BEGIN_RCPP
	NumericVector longpars(Rlongpars);
    List pars(Rpars);
    NumericVector ngroups(Rngroups);
    NumericVector J(RJ);
    int ind = 0;

    for(int g = 0; g < ngroups[0]; ++g){
        List glist = pars[g];
        for(int i = 0; i < (J[0]+1); ++i){
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
    
    IntegerMatrix fulldata(Rfulldata);
    NumericMatrix itemtrace0(Ritemtrace0);    
    NumericMatrix itemtrace1(Ritemtrace1);    
    NumericVector log_den0(Rlog_den0);
    NumericVector log_den1(Rlog_den1);
    List ret;
    NumericVector Sum0(fulldata.nrow()), Sum1(fulldata.nrow());;
    
    
    for(int i = 0; i < fulldata.nrow(); ++i){
        double rs0 = 0.0;
        double rs1 = 0.0;
        for(int j = 0; j < fulldata.ncol(); ++j){
            if(fulldata(i,j)){
                rs0 += log(itemtrace0(i,j));
                rs1 += log(itemtrace1(i,j));
            }
        }
        Sum0(i) = rs0 + log_den0(i);
        Sum1(i) = rs1 + log_den1(i);
    }
	
    ret["total_0"] = Sum0;
    ret["total_1"] = Sum1;
    return(ret);
	END_RCPP
}

