// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "Rfast.h"

using namespace arma;
using namespace Rcpp;

#ifndef MN
#define MN

using std::string;
using std::vector;

double sum_max_elems(colvec, colvec);
colvec max_elems(colvec, colvec);
mat colMaxElems(mat, colvec);

template<class T>
T colSumMaxs(mat &x, colvec y){
	T res(x.n_cols);
	for(unsigned int i=0;i<x.n_cols;++i){
		res[i] = sum_max_elems(x.col(i),y);
	}
	return res;
}

colvec euclidean_norm(mat &);
rowvec operator/(colvec x, double s);
bool my_compare_order_second(const pr<double, int> &, const pr<double, int> &);
NumericMatrix design_matrix_regr(CharacterVector x);
vec regression_only(mat, colvec);
double regression_only_col(colvec, colvec &);
double digamma(double);
double trigamma(double);
void i4mat_floyd(int, NumericVector &);
void i4mat_floyd_with_paths(const int, NumericVector &, NumericVector &);
rowvec colMedians(mat);
void combn(arma::vec &vals, const int n, const unsigned int start_idx,
		   std::vector<double> &combn_data, double *&combn_col);
int my_round(const double);
double my_round_gen_na_rm(double, const int &);
double my_round_gen_simple(double, const int &);
int len_sort_unique_int(IntegerVector);
uvec Order_rmdp(colvec &);
rowvec colvar_rmdp(mat &);
umat design_matrix_helper_big(CharacterVector);
NumericVector minus_mean(NumericVector &, const double);
void minus_c(double f[], double &, double *, int, int &);
int True(int *, int *);
bool my_all(int *, int *);
bool my_any(int *, int *);
double total_dista(NumericMatrix, NumericMatrix, const bool);
colvec pnormc(colvec);
double sum_abs(mat, mat);
NumericVector toNumbers(string, const string);
IntegerVector combine(IntegerVector, IntegerVector);
double total_euclidean_dist(NumericMatrix, const bool);
NumericMatrix euclidean_dist(NumericMatrix, const bool);
icolvec get_k_indices(rowvec, const int &);
colvec get_k_values(rowvec, const int &);
bool check_if_is_finite(double);
IntegerVector Order(NumericVector, const bool, const bool);
NumericVector Rank(NumericVector, string, const bool, const bool);
double calcDevRes(colvec, colvec, colvec);

#endif