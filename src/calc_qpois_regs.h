#ifndef _calc_qpois_regs_h_
#define _calc_qpois_regs_h_

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace std;
using namespace arma;
using namespace Rcpp;

arma::vec calc_qpois_regs(arma::mat& x, arma::vec& y, const double tol, const double ylogy, const double my);

#endif
