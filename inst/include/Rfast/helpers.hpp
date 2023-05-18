
#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;
using namespace Rcpp;

template <class T, Mfunction<void, typename T::iterator, typename T::iterator> Function, class FF>
void setResult(List &f, const int i, FF c)
{
	T y = clone(as<T>(c));
	Function(y.begin(), y.end());
	f[i] = T(y.begin(), y.end());
}

template <class T,
		  Mfunction<
			  void,
			  typename T::iterator,
			  typename T::iterator,
			  Mfunction<
				  bool,
				  typename std::remove_reference<typename T::value_type>::type,
				  typename std::remove_reference<typename T::value_type>::type>>
			  Function,
		  class FF,
		  class F>
void setResult(List &f, const int i, FF c, F cmp)
{
	T y = clone(as<T>(c));
	Function(y.begin(), y.end(), cmp);
	f[i] = y;
}

template <class T, class HELPER,
		  Mfunction<
			  void,
			  typename T::iterator,
			  typename T::iterator,
			  Mfunction<
				  bool,
				  typename std::remove_reference<typename T::value_type>::type,
				  typename std::remove_reference<typename T::value_type>::type>>
			  Function,
		  class FF,
		  class F>
void setResultParallelSection(List &f, FF s, F cmp)
{
	T y;
	int i;
#pragma omp critical
	{
		HELPER yy = *s;
		y = as<T>(yy);
		i = s - f.begin();
	}
	Function(y.begin(), y.end(), cmp);
#pragma omp critical
	{
		f[i] = HELPER(y.begin(), y.end());
	}
}

template <class T, class HELPER, Mfunction<void, typename T::iterator, typename T::iterator> Function, class FF>
void setResultParallelSection(List &f, FF s)
{
	T y;
	int i;
#pragma omp critical
	{
		HELPER yy = *s;
		y = as<T>(yy);
		i = s - f.begin();
	}
	Function(y.begin(), y.end());
#pragma omp critical
	{
		f[i] = HELPER(y.begin(), y.end());
	}
}

inline long long int get_current_nanoseconds()
{
	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

#endif