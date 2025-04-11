
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "cts.h"
#include "sw_regs.h"
#include "pc_skel.h"
#include "k_nn.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::vec perm_cor(arma::vec x, arma::vec y, const unsigned int r) {
	return calc_perm_cor(x, y, r);
}

RcppExport SEXP Rfast_perm_cor(SEXP xSEXP,SEXP ySEXP,SEXP rSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< vec >::type x(xSEXP);
    traits::input_parameter< vec >::type y(ySEXP);
    traits::input_parameter< const unsigned int >::type r(rSEXP);
    __result = perm_cor(x,y,r);
    return __result;
END_RCPP
}	

// [[Rcpp::export]]
Rcpp::NumericMatrix bic_fs_reg(Rcpp::NumericVector y, Rcpp::NumericMatrix ds, const double tol, const string type) {
	return calc_bic_fs_reg(y, ds, tol, type);
}

RcppExport SEXP Rfast_bic_fs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP tolSEXP,SEXP typeSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< NumericMatrix>::type ds(dsSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    __result = bic_fs_reg(y,ds,tol,type);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
Rcpp::NumericMatrix fs_reg(Rcpp::NumericVector y, Rcpp::NumericMatrix ds, const double sig, 
		const double tol, const string type) {
	if (!type.compare("logistic") || !type.compare("poisson")) {
		return calc_fs_reg_st(y, ds, sig, tol, type);
	}
	else if (!type.compare("quasilogistic") || !type.compare("quasipoisson")) {
		return calc_fs_reg_ext(y, ds, sig, tol, type);
	}
	Rcpp::stop("Unrecognised type.\n");
}


RcppExport SEXP Rfast_fs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP sigSEXP,SEXP tolSEXP,SEXP methodSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< NumericMatrix>::type ds(dsSEXP);
    traits::input_parameter< const double >::type sig(sigSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    __result = fs_reg(y,ds,sig,tol,method);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
Rcpp::List bs_reg(arma::vec y, arma::mat ds, const double sig, const std::string type) {
	return calc_bs_reg(y, ds, sig, type);
}

RcppExport SEXP Rfast_bs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP sigSEXP,SEXP typeSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< arma::vec >::type y(ySEXP);
    traits::input_parameter< arma::mat >::type ds(dsSEXP);
    traits::input_parameter< const double >::type sig(sigSEXP);
    traits::input_parameter< const std::string >::type type(typeSEXP);
    __result = bs_reg(y,ds,sig,type);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
Rcpp::List pc_skel(arma::mat ds, const string method, const double sig, const unsigned int r, 
		arma::mat stats_init, arma::mat pvalues_init, arma::ivec is_init_vals) {
	return calc_pc_skel(ds, method, sig, r, stats_init, pvalues_init, is_init_vals);
}

RcppExport SEXP Rfast_pc_skel(SEXP dsSEXP,SEXP methodSEXP,SEXP sigSEXP,SEXP rSEXP,SEXP stats_initSEXP,SEXP pvalues_initSEXP,SEXP is_init_valsSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< mat >::type ds(dsSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const double >::type sig(sigSEXP);
    traits::input_parameter< const unsigned int >::type r(rSEXP);
    traits::input_parameter< mat >::type stats_init(stats_initSEXP);
    traits::input_parameter< mat >::type pvalues_init(pvalues_initSEXP);
    traits::input_parameter< ivec >::type is_init_vals(is_init_valsSEXP);
    __result = pc_skel(ds,method,sig,r,stats_init,pvalues_init,is_init_vals);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
arma::mat k_nn(arma::mat ds_extra, arma::vec y, arma::mat ds, arma::uvec idxs, const std::string dist_type, const std::string type, const std::string method,
		const unsigned int freq_option, const bool mem_eff) {
	idxs -= 1;
	return calc_k_nn(ds_extra, y, ds, idxs, dist_type, type, method, freq_option, mem_eff);
}

RcppExport SEXP Rfast_k_nn(SEXP ds_extraSEXP,SEXP ySEXP,SEXP dsSEXP,SEXP idxsSEXP,SEXP dist_typeSEXP,SEXP typeSEXP,SEXP methodSEXP,SEXP freq_optionSEXP,SEXP mem_effSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< mat >::type ds_extra(ds_extraSEXP);
    traits::input_parameter< vec >::type y(ySEXP);
    traits::input_parameter< mat >::type ds(dsSEXP);
    traits::input_parameter< uvec >::type idxs(idxsSEXP);
    traits::input_parameter< const string >::type dist_type(dist_typeSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const unsigned int >::type freq_option(freq_optionSEXP);
    traits::input_parameter< const bool >::type mem_eff(mem_effSEXP);
    __result = k_nn(ds_extra,y,ds,idxs,dist_type,type,method,freq_option,mem_eff);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
Rcpp::List k_nn_cv(Rcpp::List folds, arma::vec y, arma::mat ds, arma::uvec idxs, const std::string dist_type, const std::string type, const std::string method,
		const unsigned int freq_option, const bool pred_ret, const bool mem_eff) { 
	return calc_k_nn_cv(folds, y, ds, idxs, dist_type, type, method, freq_option, pred_ret, mem_eff);
}

RcppExport SEXP Rfast_k_nn_cv(SEXP foldsSEXP,SEXP ySEXP,SEXP dsSEXP,SEXP idxsSEXP,SEXP dist_typeSEXP,SEXP typeSEXP,SEXP methodSEXP,SEXP freq_optionSEXP,SEXP pred_retSEXP,SEXP mem_effSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< List >::type folds(foldsSEXP);
    traits::input_parameter< vec >::type y(ySEXP);
    traits::input_parameter< mat >::type ds(dsSEXP);
    traits::input_parameter< uvec >::type idxs(idxsSEXP);
    traits::input_parameter< const string >::type dist_type(dist_typeSEXP);
    traits::input_parameter< const string >::type type(typeSEXP);
    traits::input_parameter< const string >::type method(methodSEXP);
    traits::input_parameter< const unsigned int >::type freq_option(freq_optionSEXP);
    traits::input_parameter< const bool >::type pred_ret(pred_retSEXP);
    traits::input_parameter< const bool >::type mem_eff(mem_effSEXP);
    __result = k_nn_cv(folds,y,ds,idxs,dist_type,type,method,freq_option,pred_ret,mem_eff);
    return __result;
END_RCPP
}
