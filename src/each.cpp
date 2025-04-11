

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

#include "Rfast.h"
#include "Rfast/types.hpp"

using namespace Rfast::Type;

using namespace arma;
using namespace Rcpp;
using namespace std;

static NumericVector eachcol_med_sum(SEXP X, SEXP Y, SEXP ind) {
	NumericMatrix x(X);
	NumericVector y(Y);
	const bool is_ind_null = Rf_isNull(ind);
	const int n = is_ind_null ? x.ncol() : LENGTH(ind);
	NumericVector f(n), tmp(x.nrow());
	if (is_ind_null) {
		for (int i = 0; i < n; ++i) {
			tmp = x.column(i) + y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	} else {
		IntegerVector indd(ind);
		for (int i = 0; i < n; ++i) {
			tmp = x.column(indd[i] - 1) + y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	}
	return f;
}

static NumericVector eachcol_med_mult(SEXP X, SEXP Y, SEXP ind) {
	NumericMatrix x(X);
	NumericVector y(Y);
	const bool is_ind_null = Rf_isNull(ind);
	const int n = is_ind_null ? x.ncol() : LENGTH(ind);
	NumericVector f(n), tmp(x.nrow());
	if (is_ind_null) {
		for (int i = 0; i < n; ++i) {
			tmp = x.column(i) * y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	} else {
		IntegerVector indd(ind);
		for (int i = 0; i < n; ++i) {
			tmp = x.column(indd[i] - 1) * y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	}
	return f;
}

static NumericVector eachcol_med_min(SEXP X, SEXP Y, SEXP ind) {
	NumericMatrix x(X);
	NumericVector y(Y);
	const bool is_ind_null = Rf_isNull(ind);
	const int n = is_ind_null ? x.ncol() : LENGTH(ind);
	NumericVector f(n), tmp(x.nrow());
	if (is_ind_null) {
		for (int i = 0; i < n; ++i) {
			tmp = x.column(i) - y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	} else {
		IntegerVector indd(ind);
		for (int i = 0; i < n; ++i) {
			tmp = x.column(indd[i] - 1) - y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	}
	return f;
}

static NumericVector eachcol_med_div(SEXP X, SEXP Y, SEXP ind) {
	NumericMatrix x(X);
	NumericVector y(Y);
	const bool is_ind_null = Rf_isNull(ind);
	const int n = is_ind_null ? x.ncol() : LENGTH(ind);
	NumericVector f(n), tmp(x.nrow());
	if (is_ind_null) {
		for (int i = 0; i < n; ++i) {
			tmp = x.column(i) / y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	} else {
		IntegerVector indd(ind);
		for (int i = 0; i < n; ++i) {
			tmp = x.column(indd[i] - 1) / y;
			f[i] = med_helper<NumericVector>(tmp.begin(), tmp.end());
		}
	}
	return f;
}

/*template<class T=NumericVector,Mfunction<T,T,T> Func>
static NumericVector eachcol_med_helper(NumericMatrix& x,NumericVector& y,SEXP ind){
  const bool is_ind_null = Rf_isNull(ind);
  const int n = is_ind_null ? x.ncol() : LENGTH(ind);
  NumericVector f(n),tmp(x.nrow());
  if(is_ind_null){
	for(int i=0;i<n;++i){
	  tmp = x.column(i);
	  tmp = Func(tmp,y);
	  f[i]=med_helper<NumericVector>(tmp.begin(),tmp.end());
	}
  }else{
	IntegerVector indd(ind);
	for(int i=0;i<n;++i){
	  tmp = x.column(indd[i]-1);
	  tmp = Func(tmp,y);
	  f[i]=med_helper<NumericVector>(tmp.begin(),tmp.end());
	}
  }
  return f;
}*/

//[[Rcpp::export]]
SEXP eachcol_apply(NumericMatrix x, NumericVector y, SEXP ind = Rfast::R::Null, const char oper = '*',
				   const string method = "sum", const bool parallel = false) {
	if (method == "sum") {
		switch (oper) {
			case '*':
				return eachcol_apply_helper<double, mmult<double>, madd<double>>(x, y, ind, parallel);
			case '/':
				return eachcol_apply_helper<double, mdiv<double>, madd<double>>(x, y, ind, parallel);
			case '+':
				return eachcol_apply_helper<double, madd<double>, madd<double>>(x, y, ind, parallel);
			case '-':
				return eachcol_apply_helper<double, mdiff<double>, madd<double>>(x, y, ind, parallel);
			case '^':
				return eachcol_apply_helper<double, std::pow, madd<double>>(x, y, ind, parallel);
		}
	} else if (method == "median") {
		switch (oper) {
			case '*':
				return eachcol_med_mult(x, y, ind);
			case '/':
				return eachcol_med_div(x, y, ind);
			case '+':
				return eachcol_med_sum(x, y, ind);
			case '-':
				return eachcol_med_min(x, y, ind);
			case '^':
				stop("Unsupported type. Type must be numeric.");
		}
	} else if (method == "max") {
		switch (oper) {
			case '*':
				return eachcol_apply_helper<double, mmult<double>, mmax<double>>(x, y, ind, parallel);
			case '/':
				return eachcol_apply_helper<double, mdiv<double>, mmax<double>>(x, y, ind, parallel);
			case '+':
				return eachcol_apply_helper<double, madd<double>, mmax<double>>(x, y, ind, parallel);
			case '-':
				return eachcol_apply_helper<double, mdiff<double>, mmax<double>>(x, y, ind, parallel);
			case '^':
				return eachcol_apply_helper<double, std::pow, mmax<double>>(x, y, ind, parallel);
		}
	} else if (method == "min") {
		switch (oper) {
			case '*':
				return eachcol_apply_helper<double, mmult<double>, mmin<double>>(x, y, ind, parallel);
			case '/':
				return eachcol_apply_helper<double, mdiv<double>, mmin<double>>(x, y, ind, parallel);
			case '+':
				return eachcol_apply_helper<double, madd<double>, mmin<double>>(x, y, ind, parallel);
			case '-':
				return eachcol_apply_helper<double, mdiff<double>, mmin<double>>(x, y, ind, parallel);
			case '^':
				return eachcol_apply_helper<double, std::pow, mmin<double>>(x, y, ind, parallel);
		}
	}
	stop("Error: wrong operation type.\n");
	return {};
}

RcppExport SEXP Rfast_eachcol_apply(SEXP xSEXP, SEXP ySEXP, SEXP ind, SEXP operSEXP, SEXP methodSEXP,
									SEXP parallelSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericVector>::type y(ySEXP);
	traits::input_parameter<const char>::type oper(operSEXP);
	traits::input_parameter<const string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = eachcol_apply(x, y, ind, oper, method, parallel);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////

static SEXP eachrow(SEXP x, SEXP y, const char oper) {
	switch (oper) {
		case '*':
			switch (type<SEXP>(x)) {
				case Types::REAL:
					return eachrow_helper<double, double, mmult<double>, Rfast::R::Real>(x, y);
				case Types::INT:
					return eachrow_helper<int, int, mmult<int>, Rfast::R::Int>(x, y);
				default:
					stop("Unsupported type. Type must be numeric or integer.");
			}
		case '+':
			switch (type<SEXP>(x)) {
				case Types::REAL:
					return eachrow_helper<double, double, madd<double>, Rfast::R::Real>(x, y);
				case Types::INT:
					return eachrow_helper<int, int, madd<int>, Rfast::R::Int>(x, y);
				default:
					stop("Unsupported type. Type must be numeric or integer.");
			}
		case '/':
			switch (type<SEXP>(x)) {
				case Types::REAL:
					return eachrow_helper<double, double, mdiv<double>, Rfast::R::Real>(x, y);
				case Types::INT:
					return eachrow_helper<int, int, mdiv<int>, Rfast::R::Int>(x, y);
				default:
					stop("Unsupported type. Type must be numeric or integer.");
			}
		case '-':
			switch (type<SEXP>(x)) {
				case Types::REAL:
					return eachrow_helper<double, double, mdiff<double>, Rfast::R::Real>(x, y);
				case Types::INT:
					return eachrow_helper<int, int, mdiff<int>, Rfast::R::Int>(x, y);
				default:
					stop("Unsupported type. Type must be numeric or integer.");
			}
		case '^':
			switch (type<SEXP>(x)) {
				case Types::REAL:
					return eachrow_helper<double, double, std::pow, Rfast::R::Real>(x, y);
				default:
					stop("Unsupported type. Type must be numeric or integer.");
			}
		case '=':
			switch (type<SEXP>(x)) {
				case Types::LOGICAL:
					return eachrow_helper<double, double, mequal<double>, Rfast::R::Lgl>(x, y);
				default:
					stop("Unsupported type. Type must be logical.");
			}

		default:
			stop("The operation doesn't supported.");
	}
	return Rfast::R::Null;
}

static double apply_eachrow(SEXP x, SEXP y, const char oper, const string method) {
	if (method == "sum") {
		switch (oper) {
			case '*':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mmult<double>, madd<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mmult<int>, madd<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '+':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, madd<double>, madd<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, madd<int>, madd<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '/':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiv<double>, madd<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiv<int>, madd<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '-':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiff<double>, madd<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiff<int>, madd<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '^':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, std::pow, madd<double>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			default:
				stop("The operation doesn't supported.");
		}
	} else if (method == "min") {
		switch (oper) {
			case '*':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mmult<double>, mmin<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mmult<int>, mmin<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '+':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, madd<double>, mmin<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, madd<int>, mmin<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '/':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiv<double>, mmin<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiv<int>, mmin<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '-':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiff<double>, mmin<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiff<int>, mmin<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '^':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, std::pow, mmin<double>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			default:
				stop("The operation doesn't supported.");
		}
	} else if (method == "max") {
		switch (oper) {
			case '*':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mmult<double>, mmax<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mmult<int>, mmax<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '+':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, madd<double>, mmax<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, madd<int>, mmax<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '/':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiv<double>, mmax<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiv<int>, mmax<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '-':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, mdiff<double>, mmax<double>>(x, y);
					case Types::INT:
						return apply_eachrow_helper<int, mdiff<int>, mmax<int>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			case '^':
				switch (type<SEXP>(x)) {
					case Types::REAL:
						return apply_eachrow_helper<double, std::pow, mmax<double>>(x, y);
					default:
						stop("Unsupported type. Type must be numeric or integer.");
				}
			default:
				stop("The operation doesn't supported.");
		}
	}
	return 0.0;
}
//[[Rcpp::export]]
SEXP eachrow(SEXP x, SEXP y, const char oper, SEXP meth) {
	return Rf_isNull(meth) ? eachrow(x, y, oper) : wrap(apply_eachrow(x, y, oper, as<string>(meth)));
}

RcppExport SEXP Rfast_eachrow(SEXP x, SEXP y, SEXP operSEXP, SEXP method) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const char>::type oper(operSEXP);
	__result = eachrow(x, y, oper, method);
	return __result;
	END_RCPP
}
