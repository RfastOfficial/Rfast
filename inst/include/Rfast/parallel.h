

#ifndef PARALLEL_H
#define PARALLEL_H

//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <algorithm>

/* definition to expand macro then apply to pragma message */

// #define STRS(x) #x
// #define STR(x) STRS(x)
// #pragma message STR(__cplusplus)

#if __cplusplus >= 201603L && (!defined(__APPLE__) && !defined(__MACH__) && !defined(__clang__)) && defined(__GNUC__) && (__GNUC__ > 9)
#define _PARALLEL_
// #pragma message "Parallel is supported"
#include <execution>
#else
// #pragma message "Parallel is not supported"
#include <exception>
#endif

namespace Rfast
{
    inline constexpr bool isStdParallelSupported()
    {
#ifdef _PARALLEL_
        return true;
#else
        return false;
#endif
    }

    template <class T>
    void sort(T begin, T end, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::sort(std::execution::par, begin, end);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::sort(begin, end);
        }
    }

    template <class T, class Function>
    void sort(T begin, T end, Function cmp, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::sort(std::execution::par, begin, end, cmp);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::sort(begin, end, cmp);
        }
    }

    template <class T>
    void stable_sort(T begin, T end, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::stable_sort(std::execution::par, begin, end);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::stable_sort(begin, end);
        }
    }

    template <class T, class Function>
    void stable_sort(T begin, T end, Function cmp, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::stable_sort(std::execution::par, begin, end, cmp);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::stable_sort(begin, end, cmp);
        }
    }

    template <class T>
    void nth_element(T begin, T middle, T end, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::nth_element(std::execution::par, begin, middle, end);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::nth_element(begin, middle, end);
        }
    }

    template <class T, class Function>
    void nth_element(T begin, T middle, T end, Function cmp, const bool parallel = false)
    {
        if (parallel)
        {
#ifdef _PARALLEL_
            std::nth_element(std::execution::par, begin, middle, end, cmp);
#else
            throw std::runtime_error("The C++ parallel library isn't supported by your system. Please, don't use the parallel argument.");
#endif
        }
        else
        {
            std::nth_element(begin, middle, end, cmp);
        }
    }
}

#endif