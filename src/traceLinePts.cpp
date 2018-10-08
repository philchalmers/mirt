#include "Misc.h"
#include "traceLinePts.h"

void itemTrace(vector<double> &P, vector<double> &Pstar, const vector<double> &a, const double *d,
        const NumericMatrix &Theta, const double *g, const double *u, const NumericVector &ot)
{
    const int N = Theta.nrow();
    const int nfact = Theta.ncol();
    const int USEOT = ot.size() > 1;

    if((*u - *g) > 0){
        for (int i = 0; i < N; ++i){
            double z = *d;
            for (int j = 0; j < nfact; ++j)
                z += a[j] * Theta(i,j);
            if(USEOT) z += ot[i];
            if(z > ABS_MAX_Z) z = ABS_MAX_Z;
            else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
            Pstar[i] = 1.0 / (1.0 + exp(-z));
            P[i] = *g + (*u - *g) * Pstar[i];
        }
    }
}

void P_dich(vector<double> &P, const vector<double> &par, const NumericMatrix &Theta,
    const NumericVector &ot, const int &N, const int &nfact)
{
    const int len = par.size();
    const double utmp = par[len-1];
    const double gtmp = par[len-2];
    const double g = antilogit(&gtmp);
    const double u = antilogit(&utmp);
    const double d = par[len-3];
    const int USEOT = ot.size() > 1;

    if((u - g) > 0){
        for (int i = 0; i < N; ++i){
            double z = d;
            for (int j = 0; j < nfact; ++j)
                z += par[j] * Theta(i,j);
            if(USEOT) z += ot[i];
            if(z > ABS_MAX_Z) z = ABS_MAX_Z;
            else if(z < -ABS_MAX_Z) z = -ABS_MAX_Z;
            P[i+N] = g + (u - g) /(1.0 + exp(-z));
            P[i] = 1.0 - P[i + N];
        }
    }
}

void P_graded(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const int &N,
    const int &nfact, const int &nint, const int &itemexp, const int &israting)
{
    const int parsize = par.size();
    vector<double> a(nfact);
    for(int i = 0; i < nfact; ++i) a[i] = par[i];
    vector<double> d(nint, 0.0);
    if(israting){
        const double t = par[parsize-1];
        for(int i = nfact; i < parsize - 1; ++i)
            d[i - nfact] = par[i] + t;
    } else {
        for(int i = nfact; i < parsize; ++i)
            d[i - nfact] = par[i];
    }
    int notordered = 0;
    for(int i = 1; i < nint; ++i)
        notordered += d[i-1] <= d[i];
    if(notordered){
        int P_size = P.size();
        for(int i = 0; i < P_size; ++i)
            P[i] = 0.0;
    } else {
        const double nullzero = 0.0, nullone = 1.0;
        NumericMatrix Pk(N, nint + 2);

        for(int i = 0; i < N; ++i)
            Pk(i,0) = 1.0;
        for(int i = 0; i < nint; ++i){
            vector<double> tmp1(N), tmp2(N);
            itemTrace(tmp1, tmp2, a, &d[i], Theta, &nullzero, &nullone, ot);
            for(int j = 0; j < N; ++j)
                Pk(j,i+1) = tmp2[j];
        }
        if(itemexp){
            int which = N * (nint + 1) - 1;
            for(int i = (Pk.ncol()-2); i >= 0; --i){
                for(int j = (N-1); j >= 0; --j){
                    P[which] = Pk(j,i) - Pk(j,i+1);
                    if(P[which] < 1e-50) P[which] = 1e-50;
                    else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
                    --which;
                }
            }
        } else {
            int which = 0;
            for(int i = 0; i < Pk.ncol(); ++i){
                for(int j = 0; j < Pk.nrow(); ++j){
                    P[which] = Pk(j,i);
                    ++which;
                }
            }
        }
    }
}

void P_gpcmIRT(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const int &N,
    const int &nfact, const int &nint)
{
    const int parsize = par.size();
    const int ncat = parsize - 1;
    const double a = par[0];
    vector<double> b(ncat-1);
    for(int i = 1; i < (parsize-1); ++i)
        b[i-1] = par[i];
    const double c = par[parsize-1];
    NumericMatrix Pk(N, ncat);
    for (int i = 0; i < N; ++i){
        vector<double> z(ncat, 1.0);
        for (int j = 1; j < ncat; ++j)
            z[j] = a * (Theta(i,0) - b[j-1]) + c + z[j-1];
        double maxz = *std::max_element(z.begin(), z.end());
        vector<double> num(ncat, 0.0);
        double den = 0.0;
        for(int j = 0; j < ncat; ++j){
            z[j] = z[j] - maxz;
            num[j] = exp(z[j]);
            den += num[j];
        }
        for(int j = 0; j < ncat; ++j)
            Pk(i, j) = num[j] / den;
    }

    int which = 0;
    for(int i = 0; i < Pk.ncol(); ++i){
        for(int j = 0; j < Pk.nrow(); ++j){
            P[which] = Pk(j,i);
            if(P[which] < 1e-50) P[which] = 1e-50;
            else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
            ++which;
        }
    }
}

void P_nominal(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const int &N,
    const int &nfact, const int &ncat, const int &returnNum,
    const int &israting)
{
    vector<double> a(nfact), ak(ncat), d(ncat);
    for(int i = 0; i < nfact; ++i)
        a[i] = par[i];
    for(int i = 0; i < ncat; ++i){
        ak[i] = par[i + nfact];
        if(israting){
            if(i)
                d[i] = par[i + nfact + ncat] + par[par.size()-1];
        } else {
            d[i] = par[i + nfact + ncat];
        }
    }
    const int USEOT = ot.size() > 1;
    NumericMatrix Num(N, ncat);
    vector<double> z(ncat);
    vector<double> Den(N, 0.0);
    vector<double> innerprod(N, 0.0);

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < nfact; ++j)
            innerprod[i] += Theta(i,j) * a[j];
    if(USEOT){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < ncat; ++j)
                z[j] = ak[j] * innerprod[i] + d[j] + ot[j];
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }
        }
    } else {
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < ncat; ++j)
                z[j] = ak[j] * innerprod[i] + d[j];
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }
        }
    }
    int which = 0;
    if(returnNum){
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j);
                ++which;
            }
        }
    } else {
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j) / Den[i];
                if(P[which] < 1e-50) P[which] = 1e-50;
                else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
                ++which;
            }
        }
    }
}

void P_nominal2(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const NumericVector &ot, const int &N,
    const int &nfact, const int &ncat, const int &returnNum,
    const int &israting)
{
    vector<double> a(nfact), d(ncat);
    NumericMatrix ak(ncat, nfact);
    for(int i = 0; i < nfact; ++i)
        a[i] = par[i];
    int ind = nfact;
    for(int j = 0; j < nfact; ++j){
    	for(int i = 0; i < ncat; ++i){
        	ak(i, j) = par[ind];
        	++ind;
    	}
    }
    for(int i = 0; i < ncat; ++i)
        d[i] = par[i + nfact + ncat*nfact];
    const int USEOT = ot.size() > 1;
    NumericMatrix Num(N, ncat);
    vector<double> z(ncat);
    vector<double> Den(N, 0.0);

    if(USEOT){
        for(int i = 0; i < N; ++i){
            for(int k = 0; k < ncat; ++k)
        		z[k] = d[k] + ot[i];
        	for(int k = 0; k < ncat; ++k)
        		for(int j = 0; j < nfact; ++j)
        			z[k] += ak(k,j) * a[j] * Theta(i,j);
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }
        }
    } else {
        for(int i = 0; i < N; ++i){
        	for(int k = 0; k < ncat; ++k)
        		z[k] = d[k];
        	for(int k = 0; k < ncat; ++k)
        		for(int j = 0; j < nfact; ++j)
        			z[k] += ak(k,j) * a[j] * Theta(i,j);
            double maxz = *std::max_element(z.begin(), z.end());
            for(int j = 0; j < ncat; ++j){
                z[j] = z[j] - maxz;
                Num(i,j) = exp(z[j]);
                Den[i] += Num(i,j);
            }
        }
    }
    int which = 0;
    if(returnNum){
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j);
                ++which;
            }
        }
    } else {
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j) / Den[i];
                if(P[which] < 1e-50) P[which] = 1e-50;
                else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
                ++which;
            }
        }
    }
}

void P_nested(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const int &N, const int &nfact, const int &ncat,
    const int &correct)
{
    NumericVector dummy(1);
	const int par_size = par.size();
    vector<double> dpar(nfact+3), npar(par_size - (nfact + 3) + nfact, 1.0);
    for(int i = 0; i < nfact+3; ++i)
        dpar[i] = par[i];
    for(int i = nfact+3; i < par_size; ++i)
        npar[i - (nfact+3) + nfact] = par[i];
    vector<double> Pd(N*2), Pn(N*(ncat-1));
    P_dich(Pd, dpar, Theta, dummy, N, nfact);
    P_nominal(Pn, npar, Theta, dummy, N, nfact, ncat-1, 0, 0);
    NumericMatrix PD = vec2mat(Pd, N, 2);
    NumericMatrix PN = vec2mat(Pn, N, ncat-1);

    int k = 0, which = 0;
    for(int i = 0; i < ncat; ++i){
        if((i+1) == correct){
            for(int j = 0; j < N; ++j){
                P[which] = PD(j,1);
                ++which;
            }
            --k;
        } else {
            for(int j = 0; j < N; ++j){
                P[which] = PD(j,0) * PN(j,k);
                if(P[which] < 1e-50) P[which] = 1e-50;
                else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
                ++which;
            }
        }
        ++k;
    }
}

void P_comp(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const int &N, const int &nfact)
{
    vector<double> a(nfact), d(nfact);
    for(int j = 0; j < nfact; ++j){
        a[j] = par[j];
        d[j] = par[j+nfact];
    }
    const double gtmp = par[nfact*2];
    const double g = antilogit(&gtmp);
    for(int i = 0; i < N; ++i) P[i+N] = 1.0;
    for(int j = 0; j < nfact; ++j)
        for(int i = 0; i < N; ++i)
            P[i+N] = P[i+N] * (1.0 / (1.0 + exp(-(a[j] * Theta(i,j) + d[j]))));
    for(int i = 0; i < N; ++i){
        P[i+N] = g + (1.0 - g) * P[i+N];
        if(P[i+N] < 1e-50) P[i+N] = 1e-50;
        else if (P[i+N] > 1.0 - 1e-50) P[i+N] = 1.0 - 1e-50;
        P[i] = 1.0 - P[i+N];
    }
}

void P_lca(vector<double> &P, const vector<double> &par, const NumericMatrix &Theta,
    const NumericMatrix &item_Q, const int &N, const int &ncat, const int &nfact,
    const int &returnNum)
{
    NumericMatrix Num(N, ncat);
    vector<double> Den(N, 0.0);

    for(int i = 0; i < N; ++i){
        vector<double> z(ncat);
        int which_par = 0;
        for(int j = 1; j < ncat; ++j){
            double innerprod = 0.0;
            for(int p = 0; p < nfact; ++p)
                innerprod += par[p + which_par] * Theta(i, p) * item_Q(j,p);
            which_par += nfact;
            z[j] = innerprod;
        }
        double maxz = *std::max_element(z.begin(), z.end());
        for(int j = 0; j < ncat; ++j){
            z[j] = z[j] - maxz;
            Num(i,j) = exp(z[j]);
            Den[i] += Num(i,j);
        }
    }
    int which = 0;
    if(returnNum){
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j);
                ++which;
            }
        }
    } else {
        for(int j = 0; j < ncat; ++j){
            for(int i = 0; i < N; ++i){
                P[which] = Num(i,j) / Den[i];
                if(P[which] < 1e-50) P[which] = 1e-50;
                else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
                ++which;
            }
        }
    }
}

void P_ideal(vector<double> &P, const vector<double> &par, const NumericMatrix &Theta,
    const NumericVector &ot, const int &N, const int &nfact)
{
    const int len = par.size();
    for (int i = 0; i < N; ++i){
        double z = par[len-1];
        for (int j = 0; j < nfact; ++j)
            z += par[j] * Theta(i,j);
        double eta = -0.5 * (z*z);
        if(eta < -20.0) eta = -20.0;
        else if(eta > -1e-10) eta = -1e-10;
        double p = exp(eta);
        P[i+N] = p;
        P[i] = 1.0 - p;
    }
}

void monopoly_z(const double &theta, const vector<double> &b, const int &k,
    double &out)
{
    out = 0.0;
    for(int i = 0; i < 2*k+1; ++i)
        out += b[i] * pow(theta, i+1);
}

void monopoly_getb (const vector<double> &a, const int &k,
    vector<double> &b)
{
    for(int i = 0; i < 2*k+1; ++i)
        b[i] = a[i] / (i+1);
}

void monopoly_geta (const int &k, const double &alpha, const double &tau,
          const vector<double> &a, vector<double>  &newa)

{
  vector<double> t(3);
  t[0] = 1;
  t[1] = -2.0 * alpha;
  t[2] = pow(alpha, 2.0) + exp(tau);
  int indx = 0;
  int indx2 = 0;
  for(int i = 0; i < 2*k-1; ++i){
    for(int j = 0; j < 2*k+1; ++j){
      if(j >= indx && j< indx + 3){
        newa[j] += a[i] * t[indx2];
        ++indx2;
      }
    }
    ++indx;
    indx2 = 0;
  }
}

void monopoly_getarec (const int &k, const double &omega, const vector<double> &alpha,
    const vector<double> &tau, vector<double> &a)
{
  vector<double> olda(2*k + 1);
  olda[0] = exp(omega);
  for(int i = 1; i <= k; ++i){
    vector<double> newa(i*2 + 1);
    std::fill(newa.begin(), newa.end(), 0.0);
    monopoly_geta(i, alpha[i-1], tau[i-1], olda, newa);
    for (int j = 0; j < i*2+1; ++j)
        olda[j] = newa[j];
  }
  for(int i = 0; i < 2*k + 1; ++i){
      a[i] = olda[i];
  }
}

void P_monopoly(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const int &N,
    const int &nfact, const int &ncat, const int &k)
{
    const double omega = par[0];
    vector<double> xic(ncat);
    vector<double> alpha(k);
    vector<double> tau(k);
    for(int i = 1; i < ncat; ++i)
        xic[i] = par[i] + xic[i-1];
    for(int i = 0; i < k; ++i){
        alpha[i] = par[i*2 + ncat];
        tau[i] = par[i*2 + ncat + 1];
    }
    vector<double> a(2*k + 1);
    vector<double> b(2*k + 1);
    NumericMatrix Num(N, ncat);
    vector<double> Den(N);
    for (int i = 0; i < N; ++i){
        double zp = 0;
        monopoly_getarec(k, omega, alpha, tau, a);
        monopoly_getb(a, k, b);
        monopoly_z(Theta(i, 0), b, k, zp);
        vector<double> z(ncat);
        for(int j = 0; j < ncat; ++j)
            z[j] = xic[j] + zp * j;
        double maxz = *std::max_element(z.begin(), z.end());
        for(int j = 0; j < ncat; ++j){
            z[j] = z[j] - maxz;
            if(z[j] < -ABS_MAX_Z) z[j] = -ABS_MAX_Z;
            Num(i,j) = exp(z[j]);
            Den[i] += Num(i,j);
        }
    }
    int which = 0;
    for(int j = 0; j < ncat; ++j){
        for(int i = 0; i < N; ++i){
            P[which] = Num(i,j) / Den[i];
            ++which;
        }
    }
}

void P_ggum(vector<double> &P, const vector<double> &par,
    const NumericMatrix &Theta, const int &N,
    const int &nfact, const int &ncat)
{
    const int C = ncat - 1;
    const int D = nfact;
    const int M = 2 * C + 1;

    vector<double> dist(N);
    for (int i = 0; i < N; ++i){
        double sumdist = 0.0;
        for (int d = 0; d < D; ++d)
            sumdist += pow(par[d], 2) * pow(Theta(i,d) - par[D+d], 2);
        dist[i] = sqrt(sumdist);
    }

    NumericMatrix Num(N, ncat);
    vector<double> Den(N);
    for (int i = 0; i < N; ++i){
        double sumtau = 0;
        vector<double> z1(ncat), z2(ncat);
        for (int w = 0; w < ncat; ++w){
            if(w > 0){
                for (int d = 0; d < D; ++d)
                    sumtau += par[d] * par[w + 2*D - 1];
            }
            z1[w] = w * dist[i] + sumtau;
            z2[w] = (M - w) * dist[i] + sumtau;
        }
        double maxz1 = *std::max_element(z1.begin(), z1.end());
        double maxz2 = *std::max_element(z2.begin(), z2.end());
        double maxz = std::max(maxz1, maxz2);
        for(int j = 0; j < ncat; ++j){
            Num(i,j) = exp(z1[j] - maxz) + exp(z2[j] - maxz);
            Den[i] += Num(i,j);
        }
    }

    int which = 0;
    for(int j = 0; j < ncat; ++j){
       for(int i = 0; i < N; ++i){
            P[which] = Num(i,j) / Den[i];
            if(P[which] < 1e-50) P[which] = 1e-50;
            else if((1.0 - P[which]) < 1e-50) P[which] = 1.0 - 1e-50;
            ++which;
        }
    }
}

RcppExport SEXP monopolyTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rncat, SEXP Rk)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const int k = as<int>(Rk);
    const int ncat = as<int>(Rncat);
    const NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    vector<double> P(N*ncat);
    P_monopoly(P, par, Theta, N, nfact, ncat, k);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

    END_RCPP
}

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

RcppExport SEXP gpcmIRTTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritemexp, SEXP Rot)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericVector ot(Rot);
    const NumericMatrix Theta(RTheta);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int ncat = par.size() - nfact;
    vector<double> P(N * ncat);
    P_gpcmIRT(P, par, Theta, ot, N, 1, ncat-1);
    NumericMatrix ret = vec2mat(P, N, ncat);
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

RcppExport SEXP gpcmTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rot, SEXP Risrating, SEXP Rhas_mat,
                                 SEXP RreturnNum)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericMatrix Theta(RTheta);
    const int israting = as<int>(Risrating);
    const int has_mat = as<int>(Rhas_mat);
    const int returnNum = as<int>(RreturnNum);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    int ncat;
    if(has_mat)
        ncat = (par.size() - nfact)/(nfact + 1);
    else
    	ncat = (par.size() - nfact)/2;
    NumericVector ot(Rot);
    vector<double> P(N*ncat);
    if(has_mat) P_nominal2(P, par, Theta, ot, N, nfact, ncat, returnNum, israting);
    	else P_nominal(P, par, Theta, ot, N, nfact, ncat, returnNum, israting);
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

RcppExport SEXP ggumTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Rncat)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const NumericMatrix Theta(RTheta);
    const int ncat = as<int>(Rncat);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    vector<double> P(N*ncat);
    P_ggum(P, par, Theta, N, nfact, ncat);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

    END_RCPP
}

RcppExport SEXP lcaTraceLinePts(SEXP Rpar, SEXP RTheta, SEXP Ritem_Q, SEXP Rncat, SEXP RreturnNum)
{
    BEGIN_RCPP

    const vector<double> par = as< vector<double> >(Rpar);
    const int ncat = as<int>(Rncat);
    const NumericMatrix Theta(RTheta);
    const NumericMatrix item_Q(Ritem_Q);
    const int nfact = Theta.ncol();
    const int N = Theta.nrow();
    const int returnNum = as<int>(RreturnNum);
    vector<double> P(N*ncat);
    P_lca(P, par, Theta, item_Q, N, ncat, nfact, returnNum);
    NumericMatrix ret = vec2mat(P, N, ncat);
    return(ret);

    END_RCPP
}

void P_switch(vector<double> &P, const vector<double> &par,
    const NumericMatrix &theta, const NumericVector &ot,
    const int &N, const int &ncat, const int &nfact2,
    const int &k, const int &itemclass)
{
    // add traceline functions for items without pre-evaluated gradient/Hessian here
    switch(itemclass){
        case 1 : // example
            P_dich(P, par, theta, ot, N, nfact2);
            break;
        case 6 :
            P_gpcmIRT(P, par, theta, ot, N, nfact2, ncat);
            break;
        case 9 :
            P_ideal(P, par, theta, ot, N, nfact2);
            break;
        case 11 :
            P_ggum(P, par, theta, N, nfact2, ncat);
            break;
        case 12 :
            P_monopoly(P, par, theta, N, nfact2, ncat, k);
            break;
    }
}

void _computeItemTrace(vector<double> &itemtrace, const NumericMatrix &Theta,
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
    int has_mat = 0;
    int k = 1;
    if(itemclass == 8)
        correct = as<int>(item.slot("correctcat"));
    if(itemclass == 12)
        k = as<int>(item.slot("k"));
    NumericMatrix item_Q;
    if(itemclass == 10)
        item_Q = as<NumericMatrix>(item.slot("item.Q"));

    /*
        1 = dich
        2 = graded
        3 = gpcm
        4 = nominal
        5 = grsm
        6 = gpcmIRT
        7 = partcomp
        8 = nestlogit
        9 = custom....have to do in R for now
        10 = lca
        11 = ggum
        12 = monopoly
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
        	has_mat = as<int>(item.slot("mat"));
        	if(has_mat) P_nominal2(P, par, theta, ot, N, nfact2, ncat, 0, 0);
            	else P_nominal(P, par, theta, ot, N, nfact2, ncat, 0, 0);
            break;
        case 4 :
            P_nominal(P, par, theta, ot, N, nfact2, ncat, 0, 0);
            break;
        case 5 :
            P_graded(P, par, theta, ot, N, nfact2, ncat-1, 1, 1);
            break;
        case 6 :
            P_gpcmIRT(P, par, theta, ot, N, nfact2, ncat);
            break;
        case 7 :
            P_comp(P, par, theta, N, nfact2);
            break;
        case 8 :
            P_nested(P, par, theta, N, nfact2, ncat, correct);
            break;
        case 9 :
            P_ideal(P, par, theta, ot, N, nfact2);
            break;
        case 10 :
            P_lca(P, par, theta, item_Q, N, ncat, nfact2, 0);
            break;
        case 11 :
            P_ggum(P, par, theta, N, nfact2, ncat);
            break;
        case 12 :
            P_monopoly(P, par, Theta, N, nfact2, ncat, k);
            break;
        default :
            P_switch(P, par, theta, ot, N, ncat, nfact, k, itemclass);
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
