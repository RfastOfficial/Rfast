#ifndef _cts_h_
#define _cts_h_

#define ARMA_64BIT_WORD

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <unordered_map>
#include <RcppArmadillo.h>
#include "Rfast.h"
#include "cts_rf.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

// Alters ds
bool adj_med_NAs(arma::mat& ds);

// Alters ds
bool adj_freq_NAs(arma::mat& ds);

arma::mat calc_rank(arma::mat& ds);

arma::mat find_combn(arma::vec& vals, const int n);

arma::uvec sub_col_max_min(arma::mat& ds, const bool cont);

arma::mat calc_pt(arma::mat& ds, 
		const int df, const bool lower_tail, const bool log_p, 
		const double add);

arma::mat ext_cols(arma::mat& ds,  const unsigned int col_a, 
		const unsigned int col_b);

// Alters dst
void cp_lt(arma::mat& src, arma::mat& dst, const int val);

// Alters ds
void adj_diag(arma::mat& ds, const double val);

arma::mat cbind_tran_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat rbind_uniq(arma::mat& ds1, arma::mat& ds2,
		const bool ass1, const bool ass2);

arma::vec to_vec(arma::mat& ds);

std::vector<double> inter(arma::vec& vals1, arma::vec& vals2);

std::vector<unsigned int> index_row_eq(arma::mat& ds, std::vector<double>& vals);

arma::mat rm_rows(arma::mat& src, arma::uvec& rows);

arma::mat rm_rows_std(arma::mat& src, std::vector<unsigned int>& rows);

arma::mat rm_cols(arma::mat& src, arma::uvec& cols);

arma::mat order_col(arma::mat& ds, const unsigned int col);

arma::mat form_cmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat_std(arma::mat& ds, std::vector<unsigned int>& rows, 
		std::vector<unsigned int>& cols);

arma::mat merge_cols(arma::mat& ds, arma::uvec& idxs);

arma::mat form_ncolcmat(arma::vec& vals, arma::mat& ds);

arma::mat form_c2mat(arma::vec& vals1, arma::vec& vals2);

arma::uvec form_vec(arma::mat& ds, const unsigned int row, arma::uvec& cols);

arma::vec form_vec_wvals(arma::mat& ds, const unsigned int row, arma::uvec& cols,
		arma::vec& vals);

// Alters ds
arma::mat append_row(arma::mat& ds, const unsigned int row, arma::vec& vals);

std::vector<unsigned int> rsum_gt_zero_idxs(arma::mat& ds);

arma::vec form_cmat_vec(arma::mat& ds, arma::rowvec& vals);

arma::mat cbind_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat adj_cols(arma::mat& src, const unsigned int dst_ncols);

arma::vec cat_ci(const unsigned int x, const unsigned int y,
		arma::uvec& cs, arma::mat& ds, arma::uvec& type, const unsigned int r);

arma::vec calc_condi(const unsigned int pos1, const unsigned int pos2, arma::uvec& cs, 
		arma::mat& ds, arma::mat& cor_ds, const std::string method, const unsigned int r);

arma::vec calc_perm_cor(arma::vec& x, arma::vec& y, const unsigned int r);

#endif
