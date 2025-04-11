
// Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <math.h>
#include <zigg/header>
#include "Rfast.h"
#include "Random.h"

using namespace Rcpp;
using namespace arma;

static zigg::Ziggurat ziggurat;
Random::uniform<Random::real> rng(0, 1);

class Gamma {
	double shape, rate, d, c;

   public:
	Gamma() {}
	Gamma(double shape, double rate)
		: shape(shape), rate(1.0 / rate), d(this->shape - 1.0 / 3.0), c(1.0 / std::sqrt(9.0 * d)) {}

	inline double operator()() {
		// Marsaglia and Tsang method for rate >= 1
		while (true) {
			double x = ziggurat.rnorm(), x2 = x * x;
			double v = 1.0 + c * x;
			v = v * v * v;
			double u = rng();

			if (v > 0 && (u < 1.0 - 0.0331 * x2 * x2 || std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v)))) {
				return d * v * rate;
			}
		}
	}
};

class BetaOne {
   protected:
	double alpha;
	Gamma alpha_d;

   public:
	BetaOne() {}
	BetaOne(double alpha) : alpha(alpha), alpha_d(Gamma(alpha, 1.0)) {}

	inline double operator()() {
		double x = alpha_d();
		return x / (x + alpha_d());
	}
};

class Beta : public BetaOne {
	double beta;
	Gamma beta_d;

   public:
	Beta(double alpha, double beta) : BetaOne(alpha), beta(beta), beta_d(Gamma(this->beta, 1.0)) {}

	inline double operator()() {
		double x = alpha_d();
		return x / (x + beta_d());
	}
};

BetaOne *_beta;
Random::uniform<Random::real> *_rng;

#ifdef _OPENMP
#pragma omp threadprivate(_beta, _rng)
#endif

static mat rotation(colvec b) {
	double ab = b[b.n_elem - 1];
	colvec ca = -b * ab;
	++ca[ca.n_elem - 1];
	ca *= (1.0 / sqrt(sum(square(ca))));
	mat A = b * ca.t();
	A -= A.t();
	A = sqrt(1 - ab * ab) * A + (ab - 1) * (b * b.t() + ca * ca.t());
	A.diag() += 1;
	return A;
}

static void randn_z(mat &res, double mean = 0.0, double stddev = 1.0)
{
    for (auto &elem : res)
    {
        elem = ziggurat.rnorm();
    }
}

static colvec rvmf_h(unsigned int n, double ca, double d1, double x0, double m, double k, double b,
					 const bool parallel = false) {
	colvec w(n, fill::none);
	const double bp1 = 1.0 + b, bm1 = 1.0 - b;

	if (parallel) {
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			_beta = new BetaOne(m);
			_rng = new Random::uniform<Random::real>(0, 1);

#ifdef _OPENMP
#pragma omp for
#endif
			for (unsigned int i = 0; i < n; ++i) {
				double ta, u, z, tmp = 0;
				for (ta = -1000.0, u = 1.0; ta - ca < log(u);) {
					z = (*_beta)();
					u = (*_rng)();
					tmp = (1.0 - bp1 * z) / (1.0 - bm1 * z);
					ta = k * tmp + d1 * log1p(-x0 * tmp);
				}
				w[i] = tmp;
			}
			delete _rng;
			delete _beta;
		}
	} else {
		BetaOne beta(m);
		for (unsigned int i = 0; i < n; ++i) {
			double ta, u, z, tmp = 0;
			for (ta = -1000.0, u = 1.0; ta - ca < log(u);) {
				z = beta();
				u = rng();
				tmp = (1.0 - bp1 * z) / (1.0 - bm1 * z);
				ta = k * tmp + d1 * log1p(-x0 * tmp);
			}
			w[i] = tmp;
		}
	}
	return w;
}

// void rvmf(unsigned int n, colvec mu, double k, mat &out) {
// 	if (k > 0.0) {
// 		mu = mu / sqrt(sum(square(mu)));
// 		colvec ini(mu.n_elem);
// 		ini[mu.n_elem - 1] = 1;
// 		const unsigned int d1 = mu.n_elem - 1;
// 		mat v1 = randn<mat>(n, d1);
// 		mat v = v1.each_col() / sqrt(sum(square(v1), 1));
// 		double b = (-2.0 * k + sqrt(4.0 * k * k + d1 * d1)) / d1;
// 		double x0 = (1.0 - b) / (1.0 + b);
// 		double m = 0.5 * d1;
// 		double ca = k * x0 + (mu.n_elem - 1) * log(1.0 - x0 * x0);
// 		colvec w = rvmf_h(n, ca, d1, x0, m, k, b);
// 		mat S = join_rows(v.each_col() % sqrt(1.0 - square(w)), w);
// 		if (approx_equal(ini, mu, "absdiff", 1e-13)) {
// 			out = S;
// 		} else if (approx_equal(-ini, mu, "absdiff", 1e-13)) {
// 			out = -S;
// 		} else {
// 			mat A = rotation(ini, mu);
// 			out = S * A.t();
// 		}
// 	} else {
// 		mat x1 = randn<mat>(n, mu.n_elem, distr_param(0, 1));
// 		out = x1.each_col() / sqrt(sum(square(x1), 1));
// 	}
// }

void rvmf(unsigned int n, colvec mu, double k, mat &out, const bool parallel) {
	if (k > 0.0) {
		mu = mu / sqrt(sum(square(mu)));
		const unsigned int d1 = mu.n_elem - 1;
		mat S(n, d1 + 1, fill::none), S_h(S.begin(), n, d1, false);
		randn_z(S_h);
		double b = (-2.0 * k + sqrt(4.0 * k * k + d1 * d1)) / d1;
		double x0 = (1.0 - b) / (1.0 + b);
		double m = 0.5 * d1;
		double ca = k * x0 + (mu.n_elem - 1) * log1p(-x0 * x0);
		S.col(d1) = rvmf_h(n, ca, d1, x0, m, k, b, parallel);
		colvec tmp = (1.0 - square(S.col(d1)));
		colvec tmp2 = sum(square(S_h), 1);
		tmp = sqrt( tmp / tmp2);
        S_h.each_col() %= tmp;

		const double M = accu(abs(-mu(span(0, d1))));
		if (M + std::abs(1 - mu[mu.n_elem-1]) < 1e-13) {
			out = S;
		} else if (M + std::abs(1 + mu[mu.n_elem-1]) < 1e-13) {
			out = -S;
		} else {
			mat A = rotation(mu);
			if (parallel)
				out = Rfast::matrix_multiplication(S, A, false, true);
			else
				out = S * A.t();
		}
	} else {
		randn_z(out);
		out.each_col() /= sqrt(sum(square(out), 1));
	}
}

NumericMatrix rvmf(unsigned int n, NumericVector Mu, double k, const bool parallel) {
	colvec mu(Mu.begin(), Mu.size(), false);
	NumericMatrix Res(n, mu.n_elem);
	mat res(Res.begin(), n, mu.n_elem, false);
	rvmf(n, mu, k, res, parallel);
	Nullable<CharacterVector> names = Mu.names();
	if (names.isNotNull()) colnames(Res) = static_cast<CharacterVector>(Mu.names());
	return Res;
}

template <class T>
T rvonmises(unsigned int n, double m, double k, const bool rads) {
	colvec u(n, fill::none);
	const double pi = M_PI;
	const double pi_180 = pi / 180.0;
	const double pi_2 = 2.0 * pi;
	if (!rads) {
		// Rcout<<__LINE__<<endl;
		m = m * pi_180;
	}
	// Rcout<<__LINE__<<endl;
	colvec mu = {cos(m), sin(m)};
	// Rcout<<__LINE__<<endl;
	if (k > 0) {
		// Rcout<<__LINE__<<endl;
		mat x(n, mu.n_elem, fill::none);
		rvmf(n, mu, k, x, false);
		// Rcout<<__LINE__<<endl;
		u = (atan(x.col(1) / x.col(0)) + pi * (x.col(0) < 0));
		// Rcout<<__LINE__<<endl;
		u = u.for_each([&](double &v) { return fmod(v, pi_2); });
		// Rcout<<__LINE__<<endl;
	} else {
		// Rcout<<__LINE__<<endl;
		u = randu<colvec>(n, distr_param(0.0, pi_2));
		// Rcout<<__LINE__<<endl;
	}
	if (!rads) {
		// Rcout<<__LINE__<<endl;
		u = u * pi_180;
	}
	return u;
}

template <>
NumericVector rvonmises<NumericVector>(unsigned int n, double m, double k, const bool rads) {
	NumericVector Res(n);
	colvec res(Res.begin(), n, false);
	res = rvonmises<colvec>(n, m, k, rads);
	return Res;
}

NumericMatrix rvonmises(unsigned int n, NumericVector M, NumericVector K, const bool rads) {
	colvec m(M.begin(), M.size(), false), k(K.begin(), K.size(), false);
	NumericMatrix F(n, m.n_elem);
	mat f(F.begin(), n, m.n_elem, false);
	for (unsigned int i = 0; i < m.n_elem; ++i) {
		f.col(i) = rvonmises<colvec>(n, m[i], k[i], rads);
	}
	Nullable<CharacterVector> names = M.names();
	if (names.isNotNull()) colnames(F) = static_cast<CharacterVector>(M.names());
	return F;
}

RcppExport SEXP Rfast_rvmf(SEXP nSEXP, SEXP mSEXP, SEXP kSEXP, SEXP parallelSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<unsigned int>::type n(nSEXP);
	traits::input_parameter<NumericVector>::type m(mSEXP);
	traits::input_parameter<double>::type k(kSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = rvmf(n, m, k, parallel);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast_rvonmises(SEXP nSEXP, SEXP mSEXP, SEXP kSEXP, SEXP radsSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<unsigned int>::type n(nSEXP);
	traits::input_parameter<const bool>::type rads(radsSEXP);

	unsigned int lenm = Rf_length(mSEXP), lenk = Rf_length(kSEXP);
	if (lenm > 1 and lenk > 1) {
		NumericVector m(mSEXP), k(kSEXP);
		__result = rvonmises(n, m, k, rads);
	} else if (lenm == 1 and lenk == 1) {
		traits::input_parameter<double>::type m(mSEXP);
		traits::input_parameter<double>::type k(kSEXP);
		__result = rvonmises<NumericVector>(n, m, k, rads);
	} else {
		throw std::runtime_error("arguments m and k must have the same length.");
	}
	return __result;
	END_RCPP
}
