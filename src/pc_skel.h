#ifndef _pc_skel_h_
#define _pc_skel_h_

#include <iostream>
#include <string>
#include <ctime>
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "cts_rf.h"
#include "cts.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

Rcpp::List calc_pc_skel(arma::mat& ds, const std::string method, const double sig, const unsigned int r, 
		arma::mat& stats_init, arma::mat& pvalues_init, arma::ivec& is_init_vals);

#endif
