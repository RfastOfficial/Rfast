#ifndef _k_nn_h_
#define _k_nn_h_

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <random>
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "cts_rf.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

arma::mat calc_k_nn(arma::mat& ds_extra, arma::vec& y, arma::mat& ds, arma::uvec& idxs,
		const std::string dist_type, const std::string type, const std::string method,
		const unsigned int freq_option, const bool mem_eff);

Rcpp::List calc_k_nn_cv(Rcpp::List& folds, arma::vec& y, arma::mat& ds, arma::uvec& idxs, 
		const std::string dist_type, const std::string type, const std::string method,
		const unsigned int freq_option, const bool pred_ret, const bool mem_eff);
#endif
