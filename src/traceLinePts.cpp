#include"Misc.h"
#include"Ps.cpp"

RcppExport SEXP traceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rot)
{
    BEGIN_RCPP

	const vector<double> par = as< vector<double> >(Rpar);
    const NumericVector ot(Rot);
    const NumericMatrix Theta(RTheta);
    const int N = Theta.nrow();
    const int nfact = Theta.ncol();
    vector<double> P(N*2);
    P_dich(P, par, Theta, ot, N, nfact);
    NumericVector ret = vec2mat(P, N, 2);
    return(ret);

	END_RCPP
}

// graded
RcppExport SEXP gradedTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritemexp, SEXP Rot, SEXP Risrating)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericVector ot(Rot);
	const NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
	const int itemexp = as<int>(Ritemexp);
    const int israting = as<int>(Risrating);
    int nint = par.size() - nfact;
    if(israting) --nint;
    int totalcat = nint + 1;
    if(!itemexp) ++totalcat;
    vector<double> P(N * totalcat);
    P_graded(P, par, Theta, ot, N, nfact, nint, itemexp, israting);
    NumericMatrix ret = vec2mat(P, N, totalcat);
    return(ret);

	END_RCPP
}

RcppExport SEXP nominalTraceLinePts(SEXP Rpar, SEXP Rncat, SEXP RTheta, SEXP RreturnNum, SEXP Rot)
{
    BEGIN_RCPP

	const vector<double> par = as< vector<double> >(Rpar);
	const int ncat = as<int>(Rncat);
	const NumericMatrix Theta(RTheta);
    const int returnNum = as<int>(RreturnNum);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    NumericVector ot(Rot);
    vector<double> P(N*ncat);
    P_nominal(P, par, Theta, ot, N, nfact, ncat, returnNum, 0);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

	END_RCPP
}

RcppExport SEXP gpcmTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Risrating)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericMatrix Theta(RTheta);
    const int israting = as<int>(Risrating);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int ncat = (par.size() - nfact)/2;
    NumericVector ot(Rot);
    vector<double> P(N*ncat);
    P_nominal(P, par, Theta, ot, N, nfact, ncat, 0, israting);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

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
    const int N = Theta.nrow();
    vector<double> P(N*ncat);
    P_nested(P, par, Theta, N, nfact, ncat, correct);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

    END_RCPP
}

RcppExport SEXP partcompTraceLinePts(SEXP Rpar, SEXP RTheta)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    vector<double> P(N*2);
    P_comp(P, par, Theta, N, nfact);
    NumericMatrix ret = vec2mat(P, N, 2);
    return(ret);

    END_RCPP
}

static void _computeItemTrace(vector<double> &itemtrace, const NumericMatrix &Theta,
    const List &pars, const NumericVector &ot, const vector<int> &itemloc, const int &which,
    const int &nfact, const int &N, const int &USEFIXED)
{
    NumericMatrix theta = Theta;
    int nfact2 = nfact;
    S4 item = pars[which];
    int ncat = as<int>(item.slot("ncat"));
    vector<double> par = as< vector<double> >(item.slot("par"));
    vector<double> P(N*ncat);
    int itemclass = as<int>(item.slot("itemclass"));
    int correct = 0;
    if(itemclass == 8)
        correct = as<int>(item.slot("correctcat"));

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
        nfact2 = nfact + itemFD.ncol();
        NumericMatrix NewTheta(Theta.nrow(), nfact2);
        for(int i = 0; i < itemFD.ncol(); ++i)
            NewTheta(_,i) = itemFD(_,i);
        for(int i = 0; i < nfact; ++i)
            NewTheta(_,i+itemFD.ncol()) = Theta(_,i);
        theta = NewTheta;
    }
    switch(itemclass){
        case 1 :
            P_dich(P, par, theta, ot, N, nfact2);
            break;
        case 2 :
            P_graded(P, par, theta, ot, N, nfact2, ncat-1, 1, 0);
            break;
        case 3 :
            P_nominal(P, par, theta, ot, N, nfact2, ncat, 0, 0);
            break;
        case 4 :
            P_nominal(P, par, theta, ot, N, nfact2, ncat, 0, 0);
            break;
        case 5 :
            P_graded(P, par, theta, ot, N, nfact2, ncat-1, 1, 1);
            break;
        case 6 :
            P_nominal(P, par, theta, ot, N, nfact2, ncat, 0, 1);
            break;
        case 7 :
            P_comp(P, par, theta, N, nfact2);
            break;
        case 8 :
            P_nested(P, par, theta, N, nfact2, ncat, correct);
            break;
        case 9 :
            break;
        default :
            Rprintf("How in the heck did you get here from a switch statement?\n");
            break;
    }
    int where = (itemloc[which]-1) * N;
    for(int i = 0; i < N*ncat; ++i)
        itemtrace[where + i] = P[i];
}

RcppExport SEXP computeItemTrace(SEXP Rpars, SEXP RTheta, SEXP Ritemloc, SEXP Roffterm)
{
    BEGIN_RCPP

    const List pars(Rpars);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix offterm(Roffterm);
    const vector<int> itemloc = as< vector<int> >(Ritemloc);
    const int J = itemloc.size() - 1;
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    vector<double> itemtrace(N * (itemloc[J]-1));
    S4 item = pars[0];
    NumericMatrix FD = item.slot("fixed.design");
    int USEFIXED = 0;
    if(FD.nrow() > 2) USEFIXED = 1;

    for(int which = 0; which < J; ++which)
        _computeItemTrace(itemtrace, Theta, pars, offterm(_, which), itemloc,
            which, nfact, N, USEFIXED);

    NumericMatrix ret = vec2mat(itemtrace, N, itemloc[J]-1);
    return(ret);

    END_RCPP
}
