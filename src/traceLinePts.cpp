#include<Rcpp.h>
#include"Misc.h"
using namespace Rcpp;

const double ABS_MAX_Z = 30;

RcppExport SEXP traceLinePts(SEXP Rpar, SEXP RTheta, SEXP RasMatrix, SEXP Rot) 
{
    BEGIN_RCPP

	NumericVector par(Rpar);
    NumericVector ot(Rot);
    IntegerVector asMatrix(RasMatrix);
    NumericMatrix Theta(RTheta);
    
    const int len = par.length();
    NumericVector a(Theta.ncol());
    const double utmp = par(len-1);
    const double gtmp = par(len-2);
    const double g = antilogit(&gtmp);
    const double u = antilogit(&utmp);
	const double d = par(len-3);
    for(int i = 0; i < Theta.ncol(); ++i)
        a(i) = par(i);    
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
    const int USEOT = ot.length() > 1;
	NumericVector P(nquad);
    NumericVector Q(nquad);
	
	NumericVector z(nquad);	
    z.fill(d);

	//compute item trace vector
	for (int j = 0; j <	nquad; ++j){
		for (int i = 0; i <	nfact; ++i)		
			z(j) += a(i) * Theta(j,i); 
	}	
    if(USEOT){
        for (int j = 0; j < nquad; ++j)
            z(j) += ot(j);
    }
	for (int i = 0; i < nquad; ++i){ 
        if(z(i) > ABS_MAX_Z) z(i) = ABS_MAX_Z;
        else if(z(i) < -ABS_MAX_Z) z(i) = -ABS_MAX_Z;
		P(i) = g + (u - g) /(1.0 + exp(-z(i)));
	}
	
    if(asMatrix(0)){
        NumericMatrix ret(nquad, 2);
        ret(_, 0) = 1.0 - P;
        ret(_, 1) = P;
        return(ret);
    } else return(P);

	END_RCPP
}

// graded
RcppExport SEXP gradedTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritemexp, SEXP Rot, SEXP Risrating) 
{
    BEGIN_RCPP

    NumericVector par(Rpar);	
    NumericVector ot(Rot);
	NumericMatrix Theta(RTheta);
	IntegerVector itemexp(Ritemexp);
    IntegerVector israting(Risrating);
    NumericVector a(Theta.ncol());
    for(int i = 0; i < Theta.ncol(); ++i)
        a(i) = par(i);
    int ncat = par.length() - Theta.ncol();
    if(israting(0)) ncat -= 1;
    NumericVector d(ncat);        
    if(israting(0)){
        const double t = par(par.length()-1);
        for(int i = Theta.ncol(); i < par.length() - 1; ++i)
            d(i - Theta.ncol()) = par(i) + t;
    } else {        
        for(int i = Theta.ncol(); i < par.length(); ++i)
            d(i - Theta.ncol()) = par(i);
    }
    const double nullzero = 0.0, nullone = 1.0;
    const int nquad = Theta.nrow();
	NumericMatrix Pk(nquad, ncat + 2);
	NumericMatrix P(nquad, ncat + 1);

	for(int i = 0; i < nquad; ++i)
        Pk(i,0) = 1.0;
    for(int i = 0; i < ncat; ++i)
        Pk(_,i+1) = itemTrace(a, &d(i), Theta, &nullzero, &nullone, ot); 
    if(itemexp(0)){
        for(int i = (Pk.ncol()-2); i >= 0; --i)
            P(_,i) = Pk(_,i) - Pk(_,i+1);
        for(int i = 0; i < P.nrow(); ++i){
            for(int j = 0; j < P.ncol(); ++j){
                if(P(i,j) < 1e-20) P(i,j) = 1e-20;
                else if((1.0 - P(i,j)) < 1e-20) P(i,j) = 1.0 - 1e-20;        
            }
        }
        return(P);
    }

    return(Pk);
	END_RCPP
}

RcppExport SEXP nominalTraceLinePts(SEXP Ra, SEXP Rak, SEXP Rd, SEXP RTheta, 
    SEXP RreturnNum, SEXP Rot) 
{
    BEGIN_RCPP

	NumericVector a(Ra);
	NumericVector ak(Rak);
	NumericVector d(Rd);	
    NumericVector ot(Rot);
	NumericMatrix Theta(RTheta);
	IntegerVector returnNum(RreturnNum);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	const int ncat = d.length();
    const int USEOT = ot.length() > 1;
    double z;

	NumericMatrix Num(nquad, ncat);
	NumericMatrix P(nquad, ncat);
	NumericVector Den(nquad);
	NumericVector innerprod(nquad);

	for(int i = 0; i < nquad; ++i)
	    for(int j = 0; j < nfact; ++j)
	        innerprod(i) += Theta(i,j) * a(j);
    if(USEOT){
        for(int i = 0; i < nquad; ++i){
            for(int j = 0; j < ncat; ++j){
                z = ak(j) * innerprod(i) + d(j) + ot(i);
                if(z > ABS_MAX_Z) z = ABS_MAX_Z;
                else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
    	        Num(i,j) = exp(z);
                Den(i) += Num(i,j);
            }        
        }
    } else {
    	for(int i = 0; i < nquad; ++i){
    	    for(int j = 0; j < ncat; ++j){
                z = ak(j) * innerprod(i) + d(j);
                if(z > ABS_MAX_Z) z = ABS_MAX_Z;
                else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
    	        Num(i,j) = exp(z);
                Den(i) += Num(i,j);
            }        
        }
    }
    if(returnNum(0)) return(Num);
	for(int i = 0; i < nquad; ++i){
	    for(int j = 0; j < ncat; ++j)
	        P(i,j) = Num(i,j) / Den(i);
    }

    return(P);
	END_RCPP
}

RcppExport SEXP gpcmTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Risrating) 
{
    BEGIN_RCPP
    
    NumericVector par(Rpar);
    NumericVector ot(Rot);
    IntegerVector israting(Risrating);
    NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    int ncat = par.length() - nfact;
    if(israting(0)) ncat -= 1;
    NumericVector a(nfact), d(ncat), ak(ncat);
    for(int i = 0; i < nfact; ++i) a(i) = par(i);
    if(israting(0)){
        const double t = par(par.length()-1);
        for(int i = nfact+1; i < par.length() - 1; ++i)
            d(i-nfact) = par(i) + t;
    } else {
        for(int i = nfact; i < par.length(); ++i)
            d(i-nfact) = par(i);
    }
    for(int i = 0; i < ak.length(); ++i) ak(i) = i;
    IntegerVector returnNum(1);
    NumericMatrix P = nominalTraceLinePts(a, ak, d, Theta, returnNum, ot);
    
    return(P);
    END_RCPP   
}

RcppExport SEXP nestlogitTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rcorrect, SEXP Rncat) 
{
    BEGIN_RCPP
    
    NumericVector par(Rpar);    
    NumericMatrix Theta(RTheta);
    IntegerVector correct(Rcorrect);
    IntegerVector ncat(Rncat);
    const int nfact = Theta.ncol();
    NumericVector dpar(nfact+3), a(nfact), d(ncat(0)-1), ak(ncat(0)-1);
    a.fill(1.0);
    for(int i = 0; i < nfact+3; ++i)
        dpar(i) = par(i);
    for(int i = 0; i < ncat(0)-1; ++i){
        ak(i) = par(i+nfact+3);
        d(i) = par(i+nfact+2+ncat(0));
    }
    NumericVector P, isfalse(1);
    NumericMatrix Pnom, traces(Theta.nrow(), ncat(0));
    P = traceLinePts(dpar, Theta, isfalse, isfalse); 
    Pnom = nominalTraceLinePts(a, ak, d, Theta, isfalse, isfalse); 
    int k = 0;
    for(int i = 0; i < traces.ncol(); ++i){
        if((i+1) == correct(0)){
            traces(_,i) = P;
            --k;
        } else {
            traces(_,i) = (1.0 - P) * Pnom(_,k);
        }
        ++k;
    }
    
    return(traces);
    END_RCPP
}

RcppExport SEXP partcompTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP RasMatrix, SEXP Rot) 
{
    BEGIN_RCPP
    
    NumericVector par(Rpar);
    NumericVector ot(Rot);
    IntegerVector asMatrix(RasMatrix);
    NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    NumericVector a(nfact), d(nfact);
    for(int j = 0; j < nfact; ++j){
        a(j) = par(j);
        d(j) = par(j+nfact);
    }
    const double gtmp = par(nfact*2);
    const double g = antilogit(&gtmp);
    NumericVector P(Theta.nrow());
    P.fill(1.0);
    
    for(int j = 0; j < nfact; ++j)
        for(int i = 0; i < Theta.nrow(); ++i)
            P(i) = P(i) * (1.0 / (1.0 + exp(-(a(j) * Theta(i,j) + d(j)))));
    for(int i = 0; i < Theta.nrow(); ++i){    
        P(i) = g + (1.0 - g) * P(i);
        if(P(i) < 1e-20) P(i) = 1e-20;
        else if (P(i) > 1.0 - 1e-20) P(i) = 1.0 - 1e-20;
    }
    if(asMatrix(0)){
        NumericMatrix ret(Theta.nrow(), 2);
        ret(_, 0) = 1.0 - P;
        ret(_, 1) = P;
        return(ret);
    } else return(P);
    END_RCPP   
}

RcppExport SEXP computeItemTrace(SEXP Rpars, SEXP RTheta, SEXP Ritemloc, SEXP Roffterm) 
{
    BEGIN_RCPP
    
    List pars(Rpars);
    NumericMatrix Theta(RTheta), offterm(Roffterm);
    IntegerVector itemloc(Ritemloc), istrue(1), isfalse(1);
    istrue.fill(1);
    const int J = itemloc.length() - 1;
    const int nfact = Theta.ncol();
    NumericMatrix itemtrace(Theta.nrow(), itemloc(J)-1);
    int where = 0;
    S4 item = pars[0];
    NumericMatrix FD = item.slot("fixed.design");
    int USEFIXED = 0;
    if(FD.nrow() > 2) USEFIXED = 1;
    
    for(int which = 0; which < J; ++which){
        S4 item = pars[which];
        IntegerVector ncat = item.slot("ncat");
        NumericVector par = item.slot("par");
        NumericVector ot = offterm(_, which);
        NumericVector a(nfact), ak(ncat(0)), d(ncat(0));
        NumericMatrix P;
        IntegerVector itemclass = item.slot("itemclass");
        IntegerVector correct;
        
        /* 
            1 = dich
            2 = graded
            3 = gpcm
            4 = nominal
            5 = grsm
            6 = rsm
            7 = partcomp
            8 = nestlogit
            9 = custom....have to do in R for now
        */
        
        if(USEFIXED){
            NumericMatrix itemFD = item.slot("fixed.design");
            NumericMatrix NewTheta(Theta.nrow(), nfact + itemFD.ncol());
            for(int i = 0; i < FD.ncol(); ++i)
                NewTheta(_,i) = itemFD(_,i);
            for(int i = 0; i < nfact; ++i)
                NewTheta(_,i+itemFD.ncol()) = Theta(_,i);
            switch(itemclass(0)){
                case 1 :
                    P = traceLinePts(par, NewTheta, istrue, ot); 
                    break;            
                case 2 :
                    P = gradedTraceLinePts(par, NewTheta, istrue, ot, isfalse); 
                    break;                
                case 3 :
                    P = gpcmTraceLinePts(par, NewTheta, ot, isfalse);
                    break;            
                case 4 :
                    for(int i = 0; i < nfact; ++i) a(i) = par(i);
                    for(int i = 0; i < ncat(0); ++i){
                        ak(i) = par(i+nfact);
                        d(i) = par(i+nfact+ncat(0)); 
                    }
                    P = nominalTraceLinePts(a, ak, d, NewTheta, isfalse, ot);
                    break;
                case 5 :
                    P = gradedTraceLinePts(par, NewTheta, istrue, ot, istrue); 
                    break;
                case 6 :
                    P = gpcmTraceLinePts(par, NewTheta, ot, istrue);
                    break;
                case 7 :
                    P = partcompTraceLinePts(par, NewTheta, istrue, ot);            
                    break;
                case 8 :
                    correct = item.slot("correctcat");
                    P = nestlogitTraceLinePts(par, NewTheta, correct, ncat);
                    break;
                case 9 :
                    continue;
                    break;
                default : 
                    Rprintf("How in the heck did you get here from a switch statement?\n");
                    break;
            }            
        } else {    
            switch(itemclass(0)){
                case 1 :
                    P = traceLinePts(par, Theta, istrue, ot); 
                    break;            
                case 2 :
                    P = gradedTraceLinePts(par, Theta, istrue, ot, isfalse); 
                    break;                
                case 3 :
                    P = gpcmTraceLinePts(par, Theta, ot, isfalse);
                    break;            
                case 4 :
                    for(int i = 0; i < nfact; ++i) a(i) = par(i);
                    for(int i = 0; i < ncat(0); ++i){
                        ak(i) = par(i+nfact);
                        d(i) = par(i+nfact+ncat(0)); 
                    }
                    P = nominalTraceLinePts(a, ak, d, Theta, isfalse, ot);
                    break;
                case 5 :
                    P = gradedTraceLinePts(par, Theta, istrue, ot, istrue); 
                    break;
                case 6 :
                    P = gpcmTraceLinePts(par, Theta, ot, istrue);
                    break;
                case 7 :
                    P = partcompTraceLinePts(par, Theta, istrue, ot);            
                    break;
                case 8 :
                    correct = item.slot("correctcat");
                    P = nestlogitTraceLinePts(par, Theta, correct, ncat);
                    break;
                case 9 :
                    continue;
                    break;
                default : 
                    Rprintf("How in the heck did you get here from a switch statement?\n");
                    break;
            }
        }
        for(int i = 0; i < P.ncol(); ++i)
            itemtrace(_, where + i) = P(_, i);
        where += P.ncol();
    }
    
    return(itemtrace);
    END_RCPP   
}
