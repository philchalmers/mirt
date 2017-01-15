#ifndef _TRACELINEPTS_H
#define _TRACELINEPTS_H

void itemTrace(vector<double> &, vector<double> &, const vector<double> &, const double *,
        const NumericMatrix &, const double *, const double *, const NumericVector &);

void P_dich(vector<double> &, const vector<double> &, const NumericMatrix &,
    const NumericVector &, const int &, const int &);

void P_graded(vector<double> &, const vector<double> &,
    const NumericMatrix &, const NumericVector &, const int &,
    const int &, const int &, const int &, const int &);

void P_gpcmIRT(vector<double> &, const vector<double> &,
    const NumericMatrix &, const NumericVector &, const int &,
    const int &, const int &);

void P_nominal(vector<double> &, const vector<double> &,
    const NumericMatrix &, const NumericVector &, const int &,
    const int &, const int &, const int &, const int &);

void P_nominal2(vector<double> &, const vector<double> &,
    const NumericMatrix &, const NumericVector &, const int &,
    const int &, const int &, const int &, const int &);

void P_nested(vector<double> &, const vector<double> &,
    const NumericMatrix &, const int &, const int &, const int &,
    const int &);

void P_comp(vector<double> &, const vector<double> &,
    const NumericMatrix &, const int &, const int &);

void P_lca(vector<double> &, const vector<double> &,
	const NumericMatrix &, const int &, const int &, const int &, const int &);

void _computeItemTrace(vector<double> &, const NumericMatrix &,
    const List &, const NumericVector &, const vector<int> &, const int &,
    const int &, const int &, const int &);

#endif
