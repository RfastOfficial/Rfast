
// Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

static mat rotation(colvec a, colvec b)
{
  double ab = dot(a, b);
  colvec ca = a - b * ab;
  ca = ca / sqrt(sum(square(ca)));
  mat A = b * ca.t();
  A = A - A.t();
  double theta = acos(ab);
  return mat(a.n_elem, a.n_elem, fill::eye) + sin(theta) * A + (cos(theta) - 1) * (b * b.t() + ca * ca.t());
}

static colvec rvmf_h(unsigned int n, double ca, double d1, double x0, double m, double k, double b)
{
  double ta, u, z, tmp = 0;
  colvec w(n, fill::none);
  for (unsigned int i = 0; i < n; ++i)
  {
    for (ta = -1000.0, u = 1.0; ta - ca < log(u);)
    {
      z = R::rbeta(m, m);
      u = R::runif(0, 1);
      tmp = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
      ta = k * tmp + d1 * log(1 - x0 * tmp);
    }
    w[i] = tmp;
  }
  return w;
}

void rvmf(unsigned int n, colvec mu, double k, mat &out)
{
  if (k > 0.0)
  {
    mu = mu / sqrt(sum(square(mu)));
    // Rcout<<__LINE__<<endl;
    colvec ini(mu.n_elem);
    // Rcout<<__LINE__<<endl;
    ini[mu.n_elem - 1] = 1;

    const unsigned int d1 = mu.n_elem - 1;
    // Rcout<<__LINE__<<endl;
    // Rcout<<n<<" "<<d1<<endl;
    mat v1 = randn<mat>(n, d1);
    // Rcout<<__LINE__<<endl;
    mat v = v1.each_col() / sqrt(sum(square(v1), 1));
    // Rcout<<__LINE__<<endl;
    double b = (-2.0 * k + sqrt(4.0 * k * k + d1 * d1)) / d1;
    double x0 = (1.0 - b) / (1.0 + b);
    double m = 0.5 * d1;
    double ca = k * x0 + (mu.n_elem - 1) * log(1.0 - x0 * x0);
    // Rcout<<__LINE__<<endl;
    colvec w = rvmf_h(n, ca, d1, x0, m, k, b);
    // Rcout<<__LINE__<<endl;
    mat S = join_rows(v.each_col() % sqrt(1.0 - square(w)), w);
    // Rcout<<__LINE__<<endl;
    if (approx_equal(ini, mu, "absdiff", 1e-13))
    {
      // Rcout<<__LINE__<<endl;
      out = S;
    }
    else if (approx_equal(-ini, mu, "absdiff", 1e-13))
    {
      // Rcout<<__LINE__<<endl;
      out = -S;
    }
    else
    {
      // Rcout<<__LINE__<<endl;
      mat A = rotation(ini, mu);
      // Rcout<<__LINE__<<endl;
      out = S * A.t();
      // Rcout<<__LINE__<<endl;
    }
  }
  else
  {
    // Rcout<<__LINE__<<endl;
    mat x1 = randn<mat>(n, mu.n_elem, distr_param(0, 1));
    out = x1.each_col() / sqrt(sum(square(x1), 1));
    // Rcout<<__LINE__<<endl;
  }
  // Rcout<<__LINE__<<endl;
}

NumericMatrix rvmf(unsigned int n, NumericVector Mu, double k)
{
    colvec mu(Mu.begin(),Mu.size(),false);
    NumericMatrix Res(n, mu.n_elem);
    mat res(Res.begin(), n, mu.n_elem, false);
    rvmf(n,mu,k,res);
    Nullable<CharacterVector> names = Mu.names();
    if (names.isNotNull())
        colnames(Res) = static_cast<CharacterVector>(Mu.names());
    return Res;
}

template<class T>
T rvonmises(unsigned int n, double m, double k, const bool rads)
{
  colvec u(n, fill::none);
  const double pi = M_PI;
  const double pi_180 = pi / 180.0;
  const double pi_2 = 2.0 * pi;
  if (!rads)
  {
    // Rcout<<__LINE__<<endl;
    m = m * pi_180;
  }
  // Rcout<<__LINE__<<endl;
  colvec mu = {cos(m), sin(m)};
  // Rcout<<__LINE__<<endl;
  if (k > 0)
  {
    // Rcout<<__LINE__<<endl;
    mat x(n,mu.n_elem,fill::none);
    rvmf(n, mu, k, x);
    // Rcout<<__LINE__<<endl;
    u = (atan(x.col(1) / x.col(0)) + pi * (x.col(0) < 0));
    // Rcout<<__LINE__<<endl;
    u = u.for_each([&](double &v)
                   { return fmod(v, pi_2); });
    // Rcout<<__LINE__<<endl;
  }
  else
  {
    // Rcout<<__LINE__<<endl;
    u = randu<colvec>(n, distr_param(0.0, pi_2));
    // Rcout<<__LINE__<<endl;
  }
  if (!rads)
  {
    // Rcout<<__LINE__<<endl;
    u = u * pi_180;
  }
  return u;
}

template<>
NumericVector rvonmises<NumericVector>(unsigned int n, double m, double k, const bool rads){
  NumericVector Res(n);
  colvec res(Res.begin(),n,false);
  res = rvonmises<colvec>(n,m,k,rads);
  return Res;
}

NumericMatrix rvonmises(unsigned int n, NumericVector M, NumericVector K, const bool rads)
{
  colvec m(M.begin(), M.size(), false), k(K.begin(), K.size(), false);
  NumericMatrix F(n, m.n_elem);
  mat f(F.begin(), n, m.n_elem, false);
  for (unsigned int i = 0; i < m.n_elem; ++i)
  {
    f.col(i) = rvonmises<colvec>(n, m[i], k[i], rads);
  }
  Nullable<CharacterVector> names = M.names();
  if (names.isNotNull())
    colnames(F) = static_cast<CharacterVector>(M.names());
  return F;
}

RcppExport SEXP Rfast_rvmf(SEXP nSEXP, SEXP mSEXP, SEXP kSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<unsigned int>::type n(nSEXP);
  traits::input_parameter<NumericVector>::type m(mSEXP);
  traits::input_parameter<double>::type k(kSEXP);
  __result = rvmf(n, m, k);
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast_rvonmises(SEXP nSEXP, SEXP mSEXP, SEXP kSEXP, SEXP radsSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<unsigned int>::type n(nSEXP);
  traits::input_parameter<const bool>::type rads(radsSEXP);

  unsigned int lenm = Rf_length(mSEXP), lenk = Rf_length(kSEXP);
  if (lenm > 1 and lenk > 1)
  {
    NumericVector m(mSEXP), k(kSEXP);
    __result = rvonmises(n, m, k, rads);
  }
  else if (lenm == 1 and lenk == 1)
  {
    traits::input_parameter<double>::type m(mSEXP);
    traits::input_parameter<double>::type k(kSEXP);
    __result = rvonmises<NumericVector>(n, m, k, rads);
  }
  else
  {
    throw std::runtime_error("arguments m and k must have the same length.");
  }
  return __result;
  END_RCPP
}