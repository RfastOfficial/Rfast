#ifndef _calc_sw_regs_h_
#define _calc_sw_regs_h_

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "cts_rf.h"
#include "cts.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

Rcpp::NumericMatrix calc_bic_fs_reg(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double tol, const std::string type);

Rcpp::NumericMatrix calc_fs_reg_st(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double sig, const double tol, const std::string type);

Rcpp::NumericMatrix calc_fs_reg_ext(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double sig, const double tol, const std::string type);

Rcpp::List calc_bs_reg(arma::vec& y, arma::mat& ds, const double sig, const std::string type);

#endif
