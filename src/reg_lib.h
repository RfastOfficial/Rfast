//Author: Stefanos Fafalios

#ifndef _reg_lib_
#define _reg_lib_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "my_k_sorted_array.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

double most_frequent_value(vec, a_node*, const int);
double average_value(vec, a_node*, const int);
double weighted_average_value(vec, a_node*, const int);
double weighted_most_frequent_value(vec, a_node*, const int);
double vmf_mle2(double, const int, const double, const double);
colvec log1pColvec(colvec,int);
double get_geom_lik(const double,const double, const double, double *, double *,const int);
mat create_id_mat(int);
double calc_f(vec, double, vec, double, double, int);
vec gold_rat3(double, vec, vec, double, vec, const int, const double);
double calc_multinom_ini(mat,vec);
vec indexesOfNum(mat,int);
double calcSumLog(mat,vec,int);
mat colvec_mat_cbind(vec, mat);
double calcDevRes(mat,vec,mat);
List varcomps_mle(NumericVector,IntegerVector,const int,const double);
double varcomps_mle2(vec, IntegerVector,int,const double);
double rint_mle2(vec, vec, IntegerVector, const double, int);
double bic_rint_reg(mat, vec, vec, IntegerVector, const int, const double, int);
vec bic_rint_regs(mat, vec, vec, IntegerVector, const double, int, int);
List glm_poisson_2(mat, vec, const double,const double, const int);
void my_pow2(vec,double *,const double,const int);
mat varcomps_mle3(vec, IntegerVector,const int, int,const bool, const double,const int);
double spml_mle2(mat, vec, vec, vec, const int, const double, const int);
vec weibull_mle2(vec, int, const double, const int);
double calc_spml_loglik(mat::col_iterator, mat::col_iterator, vec::iterator, vec::iterator, const int);

#endif
