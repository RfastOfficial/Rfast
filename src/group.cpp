// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

#include "Rfast.h"
#include "Rfast/Set.h"
#include "Rfast/types.hpp"
#include "mn.h"

using std::string;

using namespace Rcpp;

// LogicalVector group_all(LogicalVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	pr<int, int> *y = new pr<int, int>[n];
// 	int i, c = 0, k;
// 	LogicalVector::iterator xx = x.begin();
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		k = *kk - 1;
// 		y[k].first += *xx;
// 		y[k].second++;
// 		y[k].is_good = true;
// 	}
// 	for (i = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			++c;
// 		}
// 	}
// 	LogicalVector F(c);
// 	for (i = 0, k = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			F[k++] = y[i].first == y[i].second;
// 		}
// 	}
// 	delete[] y;
// 	return F;
// }

//////////////////////////////////////////////////////////////////////////////////

// NumericMatrix group_min_max(NumericVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	const double int_max = INT_MAX;
// 	NumericVector mn(n, int_max), mx(n, (double)(INT_MIN));
// 	NumericVector::iterator xx = x.begin(), minv = mn.begin(), maxv = mx.begin();
// 	int k;
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		k = *kk - 1;
// 		mx[k] = std::max(mx[k], *xx);
// 		mn[k] = std::min(mn[k], *xx);
// 	}
// 	int count_not_zero = 0;
// 	for (; minv != mn.end(); ++minv) {
// 		if (*minv != int_max) ++count_not_zero;
// 	}
// 	NumericMatrix res(2, count_not_zero);
// 	int i = 0;
// 	for (minv = mn.begin(), maxv = mx.begin(); minv != mn.end(); ++minv, ++maxv) {
// 		if (*minv != int_max) {
// 			res(0, i) = *minv;
// 			res(1, i) = *maxv;
// 			++i;
// 		}
// 	}
// 	return res;
// }

//////////////////////////////////////////////////////////////////////////////

NumericVector group_sum(NumericVector x, IntegerVector key, SEXP minn, SEXP maxx) {
	int *mn = nullptr, *mx = nullptr;
	if (!Rf_isNull(minn)) *mn = Rf_asInteger(minn);
	if (!Rf_isNull(maxx)) *mn = Rf_asInteger(maxx);

	return group_sum_helper<NumericVector, NumericVector, IntegerVector>(x, key, mn, mx);
}

RcppExport SEXP Rfast_group_sum(SEXP xSEXP, SEXP groupSEXP, SEXP minn, SEXP maxx) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type x(xSEXP);
	traits::input_parameter<IntegerVector>::type group(groupSEXP);
	__result = group_sum(x, group, minn, maxx);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////

// NumericVector group_min(NumericVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	const double int_max = INT_MAX;
// 	NumericVector f(n, int_max);
// 	NumericVector::iterator xx = x.begin(), ff = f.begin(), rr;
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		f[*kk - 1] = std::min(f[*kk - 1], *xx);
// 	}
// 	int count_not_zero = 0;
// 	for (; ff != f.end(); ++ff) {
// 		if (*ff != int_max) ++count_not_zero;
// 	}
// 	NumericVector res(count_not_zero);
// 	for (rr = res.begin(), ff = f.begin(); ff != f.end(); ++ff) {
// 		if (*ff != int_max) *rr++ = *ff;
// 	}
// 	return res;
// }

////////////////////////////////////////////////////////////////////////

/*NumericVector group_med(NumericVector x,IntegerVector group){
  const int n=x.size(),n_1=n+1;
  IntegerVector::iterator kk=group.begin();
  NumericVector f(n);
  pr<double,int> *y=new pr<double,int>[n_1];
  pr<int,int> *ind=new pr<int,int>[n];
  int i,j=0,c=0;
  NumericVector::iterator xx=x.begin();
  for(i=0;xx!=x.end();++xx,++kk,++i){
				y[i].first=*xx;
				y[i].second=*kk-1;
  }
  y[n]=pr<double,int>();
  sort(y,y+n,my_compare_order_second);
  for(i=1;i<n_1;++i){
				if(y[j].second!=y[i].second){
				  ind[y[j].second].first=j;
				  ind[y[j].second].second=i;
				  ind[y[j].second].is_good=true;
				  ++c;
				  j=i;
				}
				f[i-1]=y[i-1].first;
  }
  NumericVector F(c);
  for(int i=0,k=0;i<n;++i){
				if(ind[i].is_good){
				  F[k++]=med_helper<NumericVector>(f.begin()+ind[i].first,f.begin()+ind[i].second);
				}
  }
  delete[] y;
  delete[] ind;
  return F;
}*/

// NumericVector group_med(NumericVector x, IntegerVector group, const int length_unique, SEXP group_maxSEXP) {
// 	int group_max;
// 	group_max = Rf_asInteger(group_maxSEXP);
// 	const int n = x.size();
// 	NumericVector f(length_unique);
// 	int i = 0;
// 	vector<vector<double>> groups(group_max, vector<double>());
// 	for (i = 0; i < n; ++i) groups[group[i] - 1].push_back(x[i]);

// 	i = 0;
// 	for (auto &gr : groups) {
// 		if (gr.size() > 0) f[i++] = med_helper<vector<double>>(gr.begin(), gr.end());
// 	}
// 	return f;
// }

///////////////////////////////////////////////////////////////////

// NumericVector group_mean(NumericVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	pr<double, int> *f = new pr<double, int>[n];
// 	NumericVector::iterator xx = x.begin(), rr;
// 	int i;
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		f[*kk - 1].first += *xx;
// 		f[*kk - 1].second++;
// 	}
// 	int count_not_zero = 0;
// 	for (i = 0; i < n; ++i) {
// 		if (f[i].second != 0) ++count_not_zero;
// 	}
// 	NumericVector res(count_not_zero);
// 	for (i = 0, rr = res.begin(); i < n; ++i) {
// 		if (f[i].second != 0) *rr++ = f[i].first / f[i].second;
// 	}
// 	delete[] f;
// 	return res;
// }

/////////////////////////////////////////////////////////////////////////////

// NumericVector group_max(NumericVector x, IntegerVector key, SEXP minn, SEXP maxx) {
// 	int mn, mx;
// 	const bool is_mn_null = Rf_isNull(minn), is_mx_null = Rf_isNull(minn);
// 	if (is_mx_null && is_mn_null) {
// 		min_max<int>(key.begin(), key.end(), mn, mx);
// 	} else if (is_mx_null) {
// 		mn = Rf_asInteger(minn);
// 		maximum<int>(key.begin(), key.end(), mx);
// 	} else if (is_mn_null) {
// 		mx = Rf_asInteger(maxx);
// 		minimum<int>(key.begin(), key.end(), mn);
// 	} else {
// 		mn = Rf_asInteger(minn);
// 		mx = Rf_asInteger(maxx);
// 	}
// 	IntegerVector::iterator kk = key.begin();
// 	NumericVector f(mx - mn + 1, double(INT_MIN));
// 	NumericVector::iterator xx = x.begin(), ff = f.begin(), rr;
// 	vector<bool>::iterator ok;
// 	int index;
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		index = *kk - mn;
// 		f[index] = std::max(f[index], *xx);
// 	}
// 	int number_of_values = 0;
// 	for (auto v : f) {
// 		if (v != INT_MIN) ++number_of_values;
// 	}
// 	NumericVector res(number_of_values);
// 	for (rr = res.begin(), ff = f.begin(); ff != f.end(); ++ff) {
// 		if (*ff != INT_MIN) *rr++ = *ff;
// 	}
// 	return res;
// }

////////////////////////////////////////////////////////////////////////////

// static double mean_ad(NumericVector::iterator first, NumericVector::iterator last) {
// 	const int n = first - last + 1;
// 	const double mn = std::accumulate(first, last, 0.0) / n;
// 	double s = 0;
// 	NumericVector::iterator ff = first;
// 	for (; ff != last; ++ff) s += abs(*ff - mn);
// 	return s / n;
// }

// static double med_ad(NumericVector::iterator first, NumericVector::iterator last) {
// 	const double md = med_helper<NumericVector>(first, last);
// 	const double center = 1.4826;
// 	NumericVector::iterator ff = first;
// 	for (; ff != last; ++ff) *ff = abs(*ff - md);
// 	return med_helper<NumericVector>(first, last) * center;
// }

// NumericVector group_mad(NumericVector x, IntegerVector group, const string method) {
// 	const int n = x.size(), n_1 = n + 1;
// 	IntegerVector::iterator kk = group.begin();
// 	NumericVector f(n);
// 	pr<double, int> *y = new pr<double, int>[n_1];
// 	pr<int, int> *ind = new pr<int, int>[n];
// 	int i, j = 0, c = 0;
// 	NumericVector::iterator xx = x.begin();
// 	for (i = 0; xx != x.end(); ++xx, ++kk, ++i) {
// 		y[i].first = *xx;
// 		y[i].second = *kk - 1;
// 	}
// 	y[n] = pr<double, int>();
// 	sort(y, y + n, my_compare_order_second);
// 	for (i = 1; i < n_1; ++i) {
// 		if (y[j].second != y[i].second) {
// 			ind[y[j].second].first = j;
// 			ind[y[j].second].second = i;
// 			ind[y[j].second].is_good = true;
// 			++c;
// 			j = i;
// 		}
// 		f[i - 1] = y[i - 1].first;
// 	}
// 	NumericVector F(c);
// 	if (method == "median") {
// 		for (int i = 0, k = 0; i < n; ++i) {
// 			if (ind[i].is_good) {
// 				F[k++] = med_ad(f.begin() + ind[i].first, f.begin() + ind[i].second);
// 			}
// 		}
// 	} else if (method == "mean") {
// 		for (int i = 0, k = 0; i < n; ++i) {
// 			if (ind[i].is_good) {
// 				F[k++] = mean_ad(f.begin() + ind[i].first, f.begin() + ind[i].second);
// 			}
// 		}
// 	}
// 	delete[] y;
// 	delete[] ind;
// 	return F;
// }

////////////////////////////////////////////////////////////////////////////

// LogicalVector group_any(LogicalVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	pr<int, int> *y = new pr<int, int>[n];
// 	int i, c = 0, k;
// 	LogicalVector::iterator xx = x.begin();
// 	for (; xx != x.end(); ++xx, ++kk) {
// 		k = *kk - 1;
// 		y[k].first += *xx;
// 		y[k].is_good = true;
// 	}
// 	for (i = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			++c;
// 		}
// 	}
// 	LogicalVector F(c);
// 	for (i = 0, k = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			F[k++] = y[i].first > 0;
// 		}
// 	}
// 	delete[] y;
// 	return F;
// }

//////////////////////////////////////////////////////////////////////////

// struct var_h {
// 	double x2;
// 	double x;
// 	int n;
// 	bool is_good;
// 	var_h() : x2(0), x(0), n(0), is_good(false) {}
// };

// NumericVector group_var(NumericVector x, IntegerVector group, SEXP maxSEXP) {
// 	int n;
// 	if (Rf_isNull(maxSEXP))
// 		maximum<int>(group.begin(), group.end(), n);
// 	else
// 		n = Rf_asInteger(maxSEXP);
// 	IntegerVector::iterator kk = group.begin();
// 	var_h *y = new var_h[n];
// 	int i, c = 0, k;
// 	double xxx;
// 	NumericVector::iterator xx = x.begin();
// 	for (i = 0; xx != x.end(); ++xx, ++kk, ++i) {
// 		k = *kk - 1;
// 		xxx = *xx;
// 		y[k].x2 += xxx * xxx;
// 		y[k].x += xxx;
// 		y[k].n++;
// 		y[k].is_good = true;
// 	}
// 	for (i = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			y[i].x *= y[i].x;
// 			++c;
// 		}
// 	}
// 	NumericVector F(c);
// 	double x2, x1;
// 	int ni;
// 	for (i = 0, k = 0; i < n; ++i) {
// 		if (y[i].is_good) {
// 			x2 = y[i].x2;
// 			x1 = y[i].x;
// 			ni = y[i].n;
// 			F[k++] = (x2 - x1 / ni) / (ni - 1);
// 		}
// 	}
// 	delete[] y;
// 	return F;
// }

#define UTYPEOF(x) ((unsigned)TYPEOF(x))

template <class T, class I, auto func>
void group_s(SEXP x, SEXP ina, SEXP &indx, const bool sorted, T init_val) {
	Group<T, I, decltype(func)> s(x, ina, func, init_val);
	indx = PROTECT(Rf_allocVector(UTYPEOF(x), s.size()));
	s.values(indx, sorted);
	Rf_copyMostAttrib(x, indx);
	UNPROTECT(1);
}

template <class T, class I, class Func>
void group_b(SEXP x, SEXP ina, SEXP &indx, const bool sorted, T init_val, Func func) {
	GroupBucket<T, I> s(x, ina, init_val);
	indx = PROTECT(Rf_allocVector(UTYPEOF(x), s.size()));
	s.values(indx, sorted, func);
	Rf_copyMostAttrib(x, indx);
	UNPROTECT(1);
}

SEXP group2(SEXP x, SEXP ina, string method, const string method_mad = "", const bool sorted = true,
			const bool std = false) {
	const Types tina = type<SEXP>(ina), tx = type<SEXP>(x);
	SEXP indx = Rfast::R::Null;
	switch (tina) {
		case Types::INT: {
			if (tx == Types::INT) {
				if (method == "sum") {
					group_s<int, int, madd<int>>(x, ina, indx, sorted, 0);
				} else if (method == "max") {
					group_s<int, int, mmax<int>>(x, ina, indx, sorted, std::numeric_limits<int>::min());
				} else if (method == "min") {
					group_s<int, int, mmin<int>>(x, ina, indx, sorted, std::numeric_limits<int>::max());
				} else if (method == "mean") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return mean(ccc);
					});
				} else if (method == "sum") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return accu(ccc);
					});
				} else if (method == "median") {
					group_b<int, int>(x, ina, indx, sorted, 0, Rfast::median<GroupBucket<int, int>::Bucket>);
				} else if (method == "var") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return Rfast::var<icolvec>(ccc, std);
					});
				} else if (method == "mad") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return Rfast::mad<icolvec>(ccc, method_mad);
					});
				} else if (method == "any") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return any(ccc);
					});
				} else if (method == "all") {
					group_b<int, int>(x, ina, indx, sorted, 0, [&](GroupBucket<int, int>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return all(ccc);
					});
				}
			} else if (tx == Types::REAL) {
				if (method == "sum") {
					group_s<double, int, madd<double>>(x, ina, indx, sorted, 0);
				} else if (method == "max") {
					group_s<double, int, mmax<double>>(x, ina, indx, sorted, std::numeric_limits<double>::min());
				} else if (method == "min") {
					group_s<double, int, mmin<double>>(x, ina, indx, sorted, std::numeric_limits<double>::max());
				} else if (method == "mean") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return mean(ccc);
					});
				} else if (method == "sum") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return accu(ccc);
					});
				} else if (method == "median") {
					group_b<double, int>(x, ina, indx, sorted, 0, Rfast::median<GroupBucket<double, int>::Bucket>);
				} else if (method == "var") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return Rfast::var<colvec>(ccc, std);
					});
				} else if (method == "mad") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return Rfast::mad<colvec>(ccc, method_mad);
					});
				} else if (method == "any") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return any(ccc);
					});
				} else if (method == "all") {
					group_b<double, int>(x, ina, indx, sorted, 0, [&](GroupBucket<double, int>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return all(ccc);
					});
				}
			}
			break;
		}
		case Types::REAL: {
			if (tx == Types::INT) {
				if (method == "sum") {
					group_s<int, double, madd<int>>(x, ina, indx, sorted, 0);
				} else if (method == "max") {
					group_s<int, double, mmax<int>>(x, ina, indx, sorted, std::numeric_limits<int>::min());
				} else if (method == "min") {
					group_s<int, double, mmin<int>>(x, ina, indx, sorted, std::numeric_limits<int>::max());
				} else if (method == "mean") {
					group_b<int, double>(x, ina, indx, sorted, 0, [&](GroupBucket<int, double>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return mean(ccc);
					});
				} else if (method == "sum") {
					group_b<int, double>(x, ina, indx, sorted, 0, [&](GroupBucket<int, double>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return accu(ccc);
					});
				} else if (method == "median") {
					group_b<int, double>(x, ina, indx, sorted, 0, Rfast::median<GroupBucket<int, double>::Bucket>);
				} else if (method == "var") {
					group_b<int, double>(x, ina, indx, sorted, 0, [&](GroupBucket<int, double>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return Rfast::var<icolvec>(ccc, std);
					});
				} else if (method == "mad") {
					group_b<int, double>(x, ina, indx, sorted, 0, [&](GroupBucket<int, double>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return Rfast::mad<icolvec>(ccc, method_mad);
					});
				} else if (method == "any") {
					group_b<int, double>(x, ina, indx, sorted, 0, [&](GroupBucket<int, double>::Bucket &s) {
						icolvec ccc(s.data(), s.size(), false);
						return any(ccc);
					});
				}
			} else if (tx == Types::REAL) {
				if (method == "sum") {
					group_s<double, double, madd<double>>(x, ina, indx, sorted, 0);
				} else if (method == "max") {
					group_s<double, double, mmax<double>>(x, ina, indx, sorted, std::numeric_limits<double>::min());
				} else if (method == "min") {
					group_s<double, double, mmin<double>>(x, ina, indx, sorted, std::numeric_limits<double>::max());
				} else if (method == "mean") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return mean(ccc);
					});
				} else if (method == "sum") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return accu(ccc);
					});
				} else if (method == "median") {
					group_b<double, double>(x, ina, indx, sorted, 0,
											Rfast::median<GroupBucket<double, double>::Bucket>);
				} else if (method == "var") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return Rfast::var<colvec>(ccc, std);
					});
				} else if (method == "mad") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return Rfast::mad<colvec>(ccc, method_mad);
					});
				} else if (method == "any") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return any(ccc);
					});
				} else if (method == "all") {
					group_b<double, double>(x, ina, indx, sorted, 0, [&](GroupBucket<double, double>::Bucket &s) {
						colvec ccc(s.data(), s.size(), false);
						return all(ccc);
					});
				}
			}
			break;
		}
		// case Types::STRING:
		// {
		//     group_s<NumericVector,SEXP>(x, ina, indx, sorted);
		//     break;
		// }
		default:
			break;
	}
	return indx;
}

///////////////////////////////////////////////////////////////////////////
// RcppExport SEXP Rfast_group(SEXP xSEXP, SEXP groupSEXP, SEXP methodSEXP, SEXP minSEXP, SEXP maxSEXP,
// 							SEXP method_madSEXP) {
// 	BEGIN_RCPP
// 	RObject __result;
// 	RNGScope __rngScope;
// 	string method = as<string>(methodSEXP);
// 	traits::input_parameter<IntegerVector>::type group(groupSEXP);

// 	if (method == "all") {
// 		traits::input_parameter<LogicalVector>::type x(xSEXP);
// 		__result = group_all(x, group, maxSEXP);
// 	} else if (method == "min.max") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_min_max(x, group, maxSEXP);
// 	} else if (method == "sum") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_sum(x, group, minSEXP, maxSEXP);
// 	} else if (method == "min") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_min(x, group, maxSEXP);
// 	} else if (method == "med") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		traits::input_parameter<const int>::type mmin(minSEXP);
// 		__result = group_med(x, group, mmin, maxSEXP);
// 	} else if (method == "me an") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_mean(x, group, maxSEXP);
// 	} else if (method == "max") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_max(x, group, minSEXP, maxSEXP);
// 	} else if (method == "mad") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		traits::input_parameter<const string>::type method_mad(method_madSEXP);
// 		__result = group_mad(x, group, method_mad);
// 	} else if (method == "any") {
// 		traits::input_parameter<LogicalVector>::type x(xSEXP);
// 		__result = group_any(x, group, maxSEXP);
// 	} else if (method == "var") {
// 		traits::input_parameter<NumericVector>::type x(xSEXP);
// 		__result = group_var(x, group, maxSEXP);
// 	}
// 	return __result;
// 	END_RCPP
// }

RcppExport SEXP Rfast_group(SEXP xSEXP, SEXP groupSEXP, SEXP methodSEXP, SEXP method_madSEXP, SEXP stdSEXP,
							SEXP sortedSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type method_mad(method_madSEXP);
	traits::input_parameter<string>::type method(methodSEXP);
	traits::input_parameter<const bool>::type sorted(sortedSEXP);
	traits::input_parameter<const bool>::type std(stdSEXP);
	__result = group2(xSEXP, groupSEXP, method, method_mad, sorted, std);
	return __result;
	END_RCPP
}
