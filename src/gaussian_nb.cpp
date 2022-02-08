


//Author: Manos Papadakis

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
NumericMatrix gaussian_nb(NumericMatrix X,NumericMatrix M,NumericMatrix S,NumericVector Dets,NumericVector Con, const int k,const bool parallel){
	mat x(X.begin(),X.nrow(),X.ncol(),false);
	mat m(M.begin(),M.nrow(),M.ncol(),false);
	mat s(S.begin(),S.nrow(),S.ncol(),false);
	vec dets(Dets.begin(),Dets.size(),false);
	vec con(Con.begin(),Con.size(),false);
	NumericMatrix Res(k,X.nrow());
	mat res(Res.begin(),Res.nrow(),Res.ncol(),false);
	if(parallel){
		#pragma omp parallel for
		for (int i=0;i<k;++i){
			rowvec mi = m.row(i);
			rowvec si = s.row(i);
			rowvec y(x.n_rows);
			auto detsi = dets[i];
			auto coni = con[i];
			for (unsigned int j=0;j<x.n_rows;++j)
				y[j] = -sum(square(x.row(j) - mi)/si) -detsi + coni;
			res.row(i) = y;
		}
	}else{
		rowvec y(x.n_rows);
		rowvec mi(m.n_cols);
		rowvec si(s.n_cols);
		for (int i=0;i<k;++i){
			mi = m.row(i);
			si = s.row(i);
			auto detsi = dets[i];
			auto coni = con[i];
			for (unsigned int j=0;j<x.n_rows;++j)
				y[j] = -sum(square(x.row(j) - mi)/si) -detsi + coni;
			res.row(i) = y;
		}
	}
	return Res;
}

RcppExport SEXP Rfast_gaussian_nb(SEXP XSEXP,SEXP MSEXP,SEXP SSEXP,SEXP DetsSEXP,SEXP ConSEXP,SEXP kSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type X(XSEXP);
    traits::input_parameter< NumericMatrix >::type M(MSEXP);
    traits::input_parameter< NumericMatrix >::type S(SSEXP);
    traits::input_parameter< NumericVector >::type Dets(DetsSEXP);
    traits::input_parameter< NumericVector >::type Con(ConSEXP);
    traits::input_parameter< const int >::type k(kSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = gaussian_nb(X,M,S,Dets,Con,k,parallel);
    return __result;
END_RCPP
}

