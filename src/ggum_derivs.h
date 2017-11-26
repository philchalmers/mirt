#ifndef _GGUM_DERIVS_H
#define _GGUM_DERIVS_H

NumericVector grad_ggum(arma::colvec, arma::mat, int, int, arma::mat);

arma::mat hess_ggum (arma::colvec, arma::mat, int, int, arma::mat);

#endif
