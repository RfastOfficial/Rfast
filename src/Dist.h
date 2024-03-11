

#ifndef DIST_H
#define DIST_H

#include <Rcpp.h>
#include <string>

using Rcpp::NumericMatrix;
using std::string;

namespace Dist
{
    template<bool sqr>
    inline double euclidean(colvec &x, colvec y)
    {
        if constexpr(sqr){
            return std::sqrt(sum(square(y - x)));
        }else{
            return sum(square(y - x));
        }
    }

    inline double manhattan(colvec &x, colvec y)
    {
        return sum(abs(x - y));
    }

    inline double chi_square(colvec &x, colvec y)
    {
        return sum(square(x - y) / (x + y));
    }
    inline double soergel(colvec &x, colvec y)
    {
        return manhattan(x, y) / sum_max_elems(x, y);
    }
    inline double kulczynski(colvec &x, colvec y)
    {
        return manhattan(x, y) / sum_min_elems(x, y);
    }
    inline double wave_hedges(colvec &x, colvec y)
    {
        return sum(abs(x - y) / max_elems(x, y));
    }
    inline double motyka(colvec &x, colvec y)
    {
        return 1.0 - sum_min_elems(x, y) / sum(x + y);
    }
    inline double harmonic_mean(colvec &x, colvec y)
    {
        return 2.0 * dot(x, y) / sum(x + y);
    }
    inline double total_variation(colvec &x, colvec y)
    {
        return 0.5 * manhattan(x, y);
    }

    inline double sorensen(colvec &x, colvec y)
    {
        return sum(abs(x - y) / (x + y));
    }

    inline double max(colvec &x, colvec y)
    {
        colvec tmp = abs(x - y);
        return tmp.at(tmp.index_max());
    }
    
    inline double min(colvec &x, colvec y)
    {
        colvec tmp = abs(x - y);
        return tmp.at(tmp.index_min());
    }
    
    inline double gower(colvec &x, colvec y, const double p)
    {
        return Dist::manhattan(x, y) * p;
    }
    
    template<bool sqr>
    inline double hellinger(colvec &x, colvec y,const double p)
    {
        if constexpr(sqr){
            return std::sqrt(sum(square(x - y)))*p;
        }else{
            return sum(square(y - x))*p;
        }
    }

}
namespace DistTotal
{
}
namespace DistVector
{
}
namespace Dista
{
}
namespace DistaIndices
{
}

NumericMatrix dist(NumericMatrix x, const string method, const bool sqr = false, const int p = 0.0);
double total_dist(NumericMatrix x, const string method, const bool sqr = false, const int p = 0.0);

double total_dista(NumericMatrix Xnew, NumericMatrix X, const string method = "",
                   const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false);
#endif