// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "mn.h"
#include "Rfast.h"
#include "Rfast/types.hpp"

using namespace Rcpp;
using namespace arma;

namespace DistaTotal
{

    double euclidean(mat &xnew, mat &x, const bool sqr, const unsigned int k)
    {
        double a = 0.0;
        if (sqr)
        {
            if (k > 0)
            {

                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += accu(get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k));
                }
            }
            else
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += sum_with<square2<double>, mat>(x.each_col() - xnew.col(i));
                }
            }
        }
        else
        {
            if (k > 0)
            {

                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += accu(get_k_values(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)), k));
                }
            }
            else
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += sum_with<std::sqrt, mat>(sum(square(x.each_col() - xnew.col(i)), 0));
                }
            }
        }
        return a;
    }

    double manhattan(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += sum_with<std::abs, mat>(x.each_col() - xnew.col(i));
            }
        }
        return a;
    }

    double sorensen(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(abs(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)));
            }
        }
        return a;
    }

    double chi_square(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(sum(square(x.each_col() - xnew.col(i)) / (x.each_col() + xnew.col(i)), 0));
            }
        }
        return a;
    }

    double cosine(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        colvec norm_xnew = euclidean_norm(xnew), norm_x = euclidean_norm(x);
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(x.each_col() % xnew.col(i), 0).t() / (norm_x * norm_xnew[i]), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(sum(x.each_col() % xnew.col(i), 0).t() / (norm_x * norm_xnew[i]));
            }
        }
        return a;
    }

    double soergel(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMaxs<rowvec>(x, xnew.col(i)), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(sum(abs(x.each_col() - xnew.col(i)), 0).t() / colSumMaxs<colvec>(x, xnew.col(i)));
            }
        }
        return a;
    }

    double kulczynski(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0) / colSumMins<rowvec>(x, xnew.col(i)), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(sum(abs(x.each_col() - xnew.col(i)), 0).t() / colSumMins<colvec>(x, xnew.col(i)));
            }
        }
        return a;
    }

    double motyka(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(1.0 - colSumMins<rowvec>(x, xnew.col(i)) / sum(abs(x.each_col() + xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(1.0 - colSumMins<colvec>(x, xnew.col(i)) / sum(abs(x.each_col() + xnew.col(i)), 0));
            }
        }
        return a;
    }

    double harmonic_mean(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(x.each_col() % xnew.col(i), 0) / sum(x.each_col() + xnew.col(i), 0), k)) * 2.0;
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(sum(x.each_col() % xnew.col(i), 0) / sum(x.each_col() + xnew.col(i), 0)) * 2.0;
            }
        }
        return a;
    }

    double hellinger(mat &xnew, mat &x, const bool sqr, const unsigned int k)
    {
        double a = 0.0;
        if (sqr)
        {
            if (k > 0)
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += accu(get_k_values(sum(square(x.each_col() - xnew.col(i)), 0), k)) * 0.5;
                }
            }
            else
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += sum_with<square2<double>, mat>(x.each_col() - xnew.col(i)) * 0.5;
                }
            }
        }
        else
        {
            const double p = 1.0 / std::sqrt(2.0);
            if (k > 0)
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += accu(get_k_values(foreach<std::sqrt, rowvec>(sum(square(x.each_col() - xnew.col(i)), 0)), k)) * p;
                }
            }
            else
            {
                for (unsigned int i = 0; i < xnew.n_cols; ++i)
                {
                    a += sum_with<std::sqrt, mat>(sum(square(x.each_col() - xnew.col(i)), 0)) * p;
                }
            }
        }
        return a;
    }

    double max(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(max(abs(x.each_col() - xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(max(abs(x.each_col() - xnew.col(i)), 0));
            }
        }
        return a;
    }

    double min(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(min(abs(x.each_col() - xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(min(abs(x.each_col() - xnew.col(i)), 0));
            }
        }
        return a;
    }

    double minkowski(mat &xnew, mat &x, const double p, const unsigned int k)
    {
        double a = 0.0;
        const double p_1 = 1.0 / p;

        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(pow(sum(pow(abs(x.each_col() - xnew.col(i)), p), 0), p_1));
            }
        }
        return a;
    }

    double canberra(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        mat x_abs = abs(x);

        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(abs(x.each_col() - xnew.col(i)) / (x_abs.each_col() + abs(xnew.col(i))));
            }
        }
        return a;
    }

    double total_variation(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k)) * 0.5;
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += sum_with<std::abs, mat>(x.each_col() - xnew.col(i)) * 0.5;
            }
        }
        return a;
    }

    double kullback_leibler(mat &xnew, mat &x, const unsigned int k, const bool parallel = false)
    {
        double a = 0.0;
        mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
        fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
        fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

        if (k > 0)
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
                double tmp = accu(get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k));
                a += tmp;
            }
        }
        else
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat m = (x.each_col() - xnew.col(i)) % (log_xx.each_col() - log_xnew.col(i));
                double tmp = sum_with_condition<double, std::isfinite>(m.begin(), m.end());
                a += tmp;
            }
        }
        return a;
    }

    double jensen_shannon(mat &xnew, mat &x, const unsigned int k, const bool parallel = false)
    {
        double a = 0.0;
        mat log_xx(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
        const double log2 = std::log(2);
        fill_with<std::log, double *, double *>(x.begin(), x.end(), log_xx.begin());
        fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());
        mat x_mod_log_xx = x % log_xx;

        if (k > 0)
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat xcolj = x.each_col() + xnew.col(i);
                mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
                mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
                double tmp = accu(get_k_values(colsum_with_condition<colvec, check_if_is_finite>(m), k));
                a += tmp;
            }
        }
        else
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat xcolj = x.each_col() + xnew.col(i);
                mat xcolj_log_xcolj = xcolj % (log2 - arma::log(xcolj));
                mat m = x_mod_log_xx + (xcolj_log_xcolj.each_col() + xnew.col(i) % log_xnew.col(i));
                double tmp = sum_with_condition<double, check_if_is_finite>(m.begin(), m.end());
                a += tmp;
            }
        }
        return a;
    }

    double bhattacharyya(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(-log(sum(sqrt(x.each_col() % xnew.col(i)), 0)), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += -sum_with<std::log, mat>(sum(sqrt(x.each_col() % xnew.col(i)), 0));
            }
        }
        return a;
    }

    double jeffries_matusita(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {

            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sqrt(2.0 - 2.0 * sum(sqrt(x.each_col() % xnew.col(i)), 0)), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += sum_with<std::sqrt, mat>(2.0 - 2.0 * sum(sqrt(x.each_col() % xnew.col(i)), 0));
            }
        }
        return a;
    }

    double itakura_saito(mat &xnew, mat &x, const unsigned int k, const bool parallel = false)
    {
        double a = 0.0;
        mat log_x(x.n_rows, x.n_cols, fill::none), log_xnew(xnew.n_rows, xnew.n_cols, fill::none);
        fill_with<std::log, double *, double *>(x.begin(), x.end(), log_x.begin());
        fill_with<std::log, double *, double *>(xnew.begin(), xnew.end(), log_xnew.begin());

        if (k > 0)
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
                double tmp = accu(get_k_values(colsum_with_condition<colvec, std::isfinite>(m), k));
                a += tmp;
            }
        }
        else
        {
#pragma omp parallel for reduction(+ : a) if (parallel)
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                mat m = (x.each_col() - xnew.col(i)) % (log_x.each_col() - log_xnew.col(i));
                double tmp = sum_with_condition<double, std::isfinite>(m.begin(), m.end());
                a += tmp;
            }
        }
        return a;
    }

    double wave_hedges(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)) / colMaxElems(x, xnew.col(i)), 0), k));
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(abs(x.each_col() - xnew.col(i)) / colMaxElems(x, xnew.col(i)));
            }
        }
        return a;
    }

    double gower(mat &xnew, mat &x, const unsigned int k)
    {
        double a = 0.0;
        const double p = 1.0 / x.n_rows;
        if (k > 0)
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(get_k_values(sum(abs(x.each_col() - xnew.col(i)), 0), k)) * p;
            }
        }
        else
        {
            for (unsigned int i = 0; i < xnew.n_cols; ++i)
            {
                a += accu(abs(x.each_col() - xnew.col(i))) * p;
            }
        }
        return a;
    }

}

//[[Rcpp::export]]
double total_dista(NumericMatrix Xnew, NumericMatrix X, const string method = "", const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false)
{
    const int n = X.ncol(), nu = Xnew.ncol();
    mat xnew(Xnew.begin(), Xnew.nrow(), nu, false), x(X.begin(), X.nrow(), n, false);
    if (method == "euclidean" || p == 2)
    {
        return DistaTotal::euclidean(xnew, x, sqr, k);
    }
    else if (method == "manhattan" || p == 1)
    {
        return DistaTotal::manhattan(xnew, x, k);
    }
    else if (method == "hellinger")
    {
        return DistaTotal::hellinger(xnew, x, sqr, k);
    }
    else if (method == "maximum")
    {
        return DistaTotal::max(xnew, x, k);
    }
    else if (method == "minimum")
    {
        return DistaTotal::min(xnew, x, k);
    }
    else if (method == "minkowski")
    {
        return DistaTotal::minkowski(xnew, x, p, k);
    }
    else if (method == "canberra")
    {
        return DistaTotal::canberra(xnew, x, k);
    }
    else if (method == "bhattacharyya")
    {
        return DistaTotal::bhattacharyya(xnew, x, k);
    }
    else if (method == "jensen_shannon")
    {
        return DistaTotal::jensen_shannon(xnew, x, k, parallel);
    }
    else if (method == "itakura_saito")
    {
        return DistaTotal::itakura_saito(xnew, x, k, parallel);
    }
    else if (method == "total_variation")
    {
        return DistaTotal::total_variation(xnew, x, k);
    }
    else if (method == "kullback_leibler")
    {
        return DistaTotal::kullback_leibler(xnew, x, k, parallel);
    }
    else if (method == "chi_square")
    {
        return DistaTotal::chi_square(xnew, x, k);
    }
    else if (method == "sorensen")
    {
        return DistaTotal::sorensen(xnew, x, k);
    }
    else if (method == "soergel")
    {
        return DistaTotal::soergel(xnew, x, k);
    }
    else if (method == "cosine")
    {
        return DistaTotal::cosine(xnew, x, k);
    }
    else if (method == "wave_hedges")
    {
        return DistaTotal::wave_hedges(xnew, x, k);
    }
    else if (method == "motyka")
    {
        return DistaTotal::motyka(xnew, x, k);
    }
    else if (method == "harmonic_mean")
    {
        return DistaTotal::harmonic_mean(xnew, x, k);
    }
    else if (method == "jeffries_matusita")
    {
        return DistaTotal::jeffries_matusita(xnew, x, k);
    }
    else if (method == "gower")
    {
        return DistaTotal::gower(xnew, x, k);
    }
    else if (method == "kulczynski")
    {
        return DistaTotal::kulczynski(xnew, x, k);
    }
    else
        stop("Unsupported Method: %s", method);
    return Rfast::NA<double>::value();
}

RcppExport SEXP Rfast_total_dista(SEXP XnewSEXP, SEXP XSEXP, SEXP methodSEXP, SEXP sqrSEXP, SEXP pSEXP, SEXP kSEXP, SEXP parallelSEXP)
{
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter<NumericMatrix>::type Xnew(XnewSEXP);
    traits::input_parameter<NumericMatrix>::type X(XSEXP);
    traits::input_parameter<const string>::type method(methodSEXP);
    traits::input_parameter<const bool>::type sqr(sqrSEXP);
    traits::input_parameter<const double>::type p(pSEXP);
    traits::input_parameter<const unsigned int>::type k(kSEXP);
    traits::input_parameter<const bool>::type parallel(parallelSEXP);

    __result = total_dista(Xnew, X, method, sqr, p, k, parallel);
    return __result;

    END_RCPP
}