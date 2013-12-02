#include"Misc.h"

RcppExport SEXP traceLinePts(SEXP Rpar, SEXP RTheta, SEXP RasMatrix, SEXP Rot) 
{
    BEGIN_RCPP

	const vector<double> par = as< vector<double> >(Rpar);
    const vector<double> ot = as< vector<double> >(Rot);
    const NumericMatrix Theta(RTheta);
    const int asMatrix = as<int>(RasMatrix);
    const int N = Theta.nrow();
    const int nfact = Theta.ncol();
    
    const int len = par.size();    
    const double utmp = par[len-1];
    const double gtmp = par[len-2];
    const double g = antilogit(&gtmp);
    const double u = antilogit(&utmp);
	const double d = par[len-3];
    const int USEOT = ot.size() > 1;
    
    vector<double> P(N);    
	for (int i = 0; i <	N; ++i){
        double z = d;        
		for (int j = 0; j <	nfact; ++j)		
			z += par[j] * Theta(i,j); 
        if(USEOT) z += ot[i];
        if(z > ABS_MAX_Z) z = ABS_MAX_Z;
        else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
        P[i] = g + (u - g) /(1.0 + exp(-z));    
	}
	
    if(asMatrix){
        NumericMatrix ret(N, 2);
        for(int j = 0; j < N; ++j){
            ret(j, 0) = 1.0 - P[j];
            ret(j, 1) = P[j];
        }
        return(ret);
    } else return(wrap(P));

	END_RCPP
}

// graded
RcppExport SEXP gradedTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritemexp, SEXP Rot, SEXP Risrating) 
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);	
    const vector<double> ot = as< vector<double> >(Rot);
	const NumericMatrix Theta(RTheta);
	const int itemexp = as<int>(Ritemexp);
    const int israting = as<int>(Risrating);
    const int parsize = par.size();
    
    vector<double> a(Theta.ncol());
    for(int i = 0; i < Theta.ncol(); ++i) a[i] = par[i];
    int ncat = parsize - Theta.ncol();
    if(israting) ncat -= 1;
    vector<double> d(ncat,0.0);        
    if(israting){
        const double t = par[parsize-1];
        for(int i = Theta.ncol(); i < parsize - 1; ++i)
            d[i - Theta.ncol()] = par[i] + t;
    } else {        
        for(int i = Theta.ncol(); i < parsize; ++i)
            d[i - Theta.ncol()] = par[i];
    }
    const double nullzero = 0.0, nullone = 1.0;
    const int nquad = Theta.nrow();
	NumericMatrix Pk(nquad, ncat + 2);
	NumericMatrix P(nquad, ncat + 1);

	for(int i = 0; i < nquad; ++i)
        Pk(i,0) = 1.0;
    for(int i = 0; i < ncat; ++i){        
        vector<double> tmp1(nquad), tmp2(nquad);
        itemTrace(tmp1, tmp2, a, &d[i], Theta, &nullzero, &nullone, ot);
        for(int j = 0; j < nquad; ++j)
            Pk(j,i+1) = tmp2[j];
    }
    if(itemexp){
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

	const vector<double> a = as< vector<double> >(Ra);
	const vector<double> ak = as< vector<double> >(Rak);
	const vector<double> d = as< vector<double> >(Rd);	
    const vector<double> ot = as< vector<double> >(Rot);
	const NumericMatrix Theta(RTheta);
	const int returnNum = as<int>(RreturnNum);
    const int nquad = Theta.nrow();
	const int nfact = Theta.ncol();
	const int ncat = d.size();
    const int USEOT = ot.size() > 1;

	NumericMatrix Num(nquad, ncat);
	NumericMatrix P(nquad, ncat);
    vector<double> z(ncat);
	vector<double> Den(nquad, 0.0);
	vector<double> innerprod(nquad, 0.0);

	for(int i = 0; i < nquad; ++i)
	    for(int j = 0; j < nfact; ++j)
	        innerprod[i] += Theta(i,j) * a[j];
    if(USEOT){
        for(int i = 0; i < nquad; ++i){
            for(int j = 0; j < ncat; ++j)
                z[j] = ak[j] * innerprod[i] + d[j] + ot[j];
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                if(z[j] < -ABS_MAX_Z) z[j] = -ABS_MAX_Z;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }       
        }
    } else {
    	for(int i = 0; i < nquad; ++i){
    	    for(int j = 0; j < ncat; ++j)
                z[j] = ak[j] * innerprod[i] + d[j];
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                if(z[j] < -ABS_MAX_Z) z[j] = -ABS_MAX_Z;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }
        }
    }
    if(returnNum) return(Num);
	for(int i = 0; i < nquad; ++i){
	    for(int j = 0; j < ncat; ++j)
	        P(i,j) = Num(i,j) / Den[i];
    }

    return(P);
	END_RCPP
}

RcppExport SEXP gpcmTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Risrating) 
{
    BEGIN_RCPP
    
    const vector<double> par = as< vector<double> >(Rpar);
    const int israting = as<int>(Risrating);
    const NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    const int parsize = par.size();
    int ncat = parsize - nfact;
    if(israting) ncat -= 1;

    NumericVector a(nfact), d(ncat), ak(ncat);
    for(int i = 0; i < nfact; ++i) a(i) = par[i];
    if(israting){
        const double t = par[parsize-1];
        for(int i = nfact+1; i < parsize - 1; ++i)
            d(i-nfact) = par[i] + t;
    } else {
        for(int i = nfact; i < parsize; ++i)
            d(i-nfact) = par[i];
    }
    for(int i = 0; i < ak.length(); ++i) ak(i) = i;
    IntegerVector returnNum(1);
    NumericMatrix P = nominalTraceLinePts(a, ak, d, Theta, returnNum, Rot);
    
    return(P);
    END_RCPP   
}

RcppExport SEXP nestlogitTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rcorrect, SEXP Rncat) 
{
    BEGIN_RCPP
    
    const vector<double> par = as< vector<double> >(Rpar);    
    const NumericMatrix Theta(RTheta);
    const int correct = as<int>(Rcorrect);
    const int ncat = as<int>(Rncat);
    const int nfact = Theta.ncol();

    NumericVector dpar(nfact+3), a(nfact, 1.0), d(ncat-1), ak(ncat-1);
    for(int i = 0; i < nfact+3; ++i)
        dpar(i) = par[i];
    for(int i = 0; i < ncat-1; ++i){
        ak(i) = par[i+nfact+3];
        d(i) = par[i+nfact+2+ncat];
    }
    const IntegerVector isfalse(1);
    NumericMatrix Pnom, traces(Theta.nrow(), ncat);
    NumericVector P = traceLinePts(dpar, Theta, isfalse, isfalse); 
    Pnom = nominalTraceLinePts(a, ak, d, Theta, isfalse, isfalse); 
    int k = 0;
    for(int i = 0; i < traces.ncol(); ++i){
        if((i+1) == correct){
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
    
    const vector<double> par = as< vector<double> >(Rpar);
    const NumericMatrix Theta(RTheta);
    const int asMatrix = as<int>(RasMatrix);
    const int nfact = Theta.ncol();
    vector<double> a(nfact), d(nfact);
    for(int j = 0; j < nfact; ++j){
        a[j] = par[j];
        d[j] = par[j+nfact];
    }
    const double gtmp = par[nfact*2];
    const double g = antilogit(&gtmp);
    vector<double> P(Theta.nrow(), 1.0);
    
    for(int j = 0; j < nfact; ++j)
        for(int i = 0; i < Theta.nrow(); ++i)
            P[i] = P[i] * (1.0 / (1.0 + exp(-(a[j] * Theta(i,j) + d[j]))));
    for(int i = 0; i < Theta.nrow(); ++i){    
        P[i] = g + (1.0 - g) * P[i];
        if(P[i] < 1e-20) P[i] = 1e-20;
        else if (P[i] > 1.0 - 1e-20) P[i] = 1.0 - 1e-20;
    }
    if(asMatrix){
        NumericMatrix ret(Theta.nrow(), 2);
        for(int j = 0; j < Theta.nrow(); ++j){
            ret(j, 0) = 1.0 - P[j];
            ret(j, 1) = P[j];
        }
        return(ret);
    } else return(wrap(P));
    
    END_RCPP   
}

RcppExport SEXP computeItemTrace(SEXP Rpars, SEXP RTheta, SEXP Ritemloc, SEXP Roffterm) 
{
    BEGIN_RCPP
    
    const List pars(Rpars);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix offterm(Roffterm);
    const IntegerVector itemloc(Ritemloc);
    const IntegerVector istrue(1, 1);
    const IntegerVector isfalse(1);
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
