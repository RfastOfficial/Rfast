

#ifndef DIST_H
#define DIST_H

#include <Rcpp.h>
#include <string>

using Rcpp::NumericMatrix;
using std::string;

namespace Dist
{
    inline double manhattan(colvec x, colvec y)
    {
        return sum(abs(x - y));
    }
    
    NumericMatrix euclidean(NumericMatrix, const bool);
}
namespace DistTotal
{
    double euclidean(NumericMatrix, const bool);
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

double total_dista(NumericMatrix Xnew, NumericMatrix X, const string method = "",
                   const bool sqr = false, const double p = 0.0, const unsigned int k = 0, const bool parallel = false);
#endif