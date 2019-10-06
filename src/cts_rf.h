#ifndef _cts_rf_h_
#define _cts_rf_h_

#include <vector>
#include <RcppArmadillo.h>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

using namespace std;
using namespace arma;
using namespace Rcpp;

double calc_med_rf(std::vector<double>& x);

arma::umat calc_dist_mem_eff_rf(arma::mat xnew, arma::mat x, const unsigned int k, const bool is_euclidean);

arma::mat calc_dist_rf(arma::mat xnew, arma::mat x, const bool is_euclidean);

Rcpp::List g2_test_univ(arma::mat& data, arma::uvec& dc);

Rcpp::List g2_test(arma::mat& data, const unsigned int x, const unsigned int y, 
		arma::uvec& cs, arma::uvec& dc);

Rcpp::List g2_test_perm(arma::mat& data, const unsigned int x, const unsigned int y,
		arma::uvec& cs, arma::uvec& dc, const unsigned int nperm);

NumericVector logistic_only(NumericMatrix& X, NumericVector& Y, const double my);

NumericVector poisson_only(NumericMatrix& X, NumericVector& Y, const double ylogy, const double my);

double glm_logistic(NumericMatrix& X, NumericVector& Y, const double my);

double arma_glm_logistic(mat x, vec y, const double my);

double glm_poisson(NumericMatrix& X, NumericVector& Y, const double ylogy, const double my);

double arma_glm_poisson(mat x, vec y, const double ylogy, const double my);

NumericVector qs_binom_only(NumericMatrix& X, NumericVector& Y, const double my);

NumericVector glm_qs_binom(NumericMatrix& X, NumericVector& Y, const double my);

NumericVector qs_poisson_only(NumericMatrix& X, NumericVector& Y, const double ylogy, const double my);

NumericVector glm_qs_poisson(NumericMatrix& X, NumericVector& Y, const double ylogy, const double my);

#endif
