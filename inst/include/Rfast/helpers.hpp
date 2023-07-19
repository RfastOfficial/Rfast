
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

inline unsigned int get_num_of_threads()
{
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 0;
#endif
}

template <class T, class HELPER, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
typename T::value_type parallelSingleIteratorWithoutCopy(FF s)
{
	T y;
#pragma omp critical
	{
		HELPER yy(s);
		y = T(yy.begin(), yy.size(), false);
	}
	return *Function(y.begin(), y.end());
}

template <class T, class Helper, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
typename T::value_type singleIteratorWithoutCopy(FF c)
{
	Helper h(c);
	T y(h.begin(), h.size(), false);
	return *Function(y.begin(), y.end());
}

template <class T, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
void singleIteratorWithoutCopy(List &f, const int i, FF c)
{
	T y(c);
	f[i] = T::create(*Function(y.begin(), y.end()));
}

template <class T, class Helper, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
void singleIteratorWithoutCopy(mat &f, const int i, FF c)
{
	Helper h(c);
	T y(h.begin(), h.size(), false);
	f[i] = *Function(y.begin(), y.end());
}

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
void setResult(mat &f, const int i, FF c, F cmp)
{
	T y = as<T>(c);
	Function(y.begin(), y.end(), cmp);
	f.col(i) = y;
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
void setResultParallelSection(mat &f, FF s, const int i, F cmp)
{
	T y;
#pragma omp critical
	{
		HELPER yy(s);
		y = as<T>(yy);
	}
	Function(y.begin(), y.end(), cmp);
	f.col(i) = y;
}

template <class T, class HELPER, Mfunction<void, typename T::iterator, typename T::iterator> Function, class FF>
void setResultParallelSection(mat &f, FF s, const int i)
{
	T y;
#pragma omp critical
	{
		HELPER yy(s);
		y = as<T>(yy);
	}
	Function(y.begin(), y.end());
	f.col(i) = y;
}

template <class T, Mfunction<void, typename T::iterator, typename T::iterator> Function, class FF>
void setResult(mat &f, const int i, FF c)
{
	T y = as<T>(c);
	Function(y.begin(), y.end());
	f.col(i) = y;
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

template <class T, class HELPER, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
void parallelSingleIteratorWithoutCopy(List &f, FF s)
{
	T y;
	int i;
#pragma omp critical
	{
		HELPER yy(s);
		y = as<T>(yy);
		i = s - f.begin();
	}
	auto value = *Function(y.begin(), y.end());
#pragma omp critical
	{
		f[i] = HELPER::create(value);
	}
}

template <class T, class HELPER, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function, class FF>
void parallelSingleIteratorWithoutCopy(colvec &f, FF s)
{
	T y;
	int i;
#pragma omp critical
	{
		HELPER yy(s);
		y = T(yy.begin(), yy.size(), false);
		i = s - f.begin();
	}
	auto value = *Function(y.begin(), y.end());
#pragma omp critical
	{
		f[i] = value;
	}
}

template <class RET, class T, class HELPER, Mfunction<std::pair<typename T::iterator, typename T::iterator>, typename T::iterator, typename T::iterator> Function, class FF>
RET parallelSingleIteratorWithoutCopy(FF s)
{
	T y;
#pragma omp critical
	{
		HELPER yy(s);
		y = T(yy.begin(), yy.size(), false);
	}
	auto v = Function(y.begin(), y.end());
	return {static_cast<typename RET::value_type>(*v.first), static_cast<typename RET::value_type>(*v.second)};
}

template <class RET, class T, class Helper, Mfunction<std::pair<typename T::iterator, typename T::iterator>, typename T::iterator, typename T::iterator> Function, class FF>
RET singleIteratorWithoutCopy(FF c)
{
	Helper h(c);
	T y(h.begin(), h.size(), false);
	auto v = Function(y.begin(), y.end());
	return {static_cast<typename RET::value_type>(*v.first), static_cast<typename RET::value_type>(*v.second)};
}

inline long long int get_current_nanoseconds()
{
	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}
#endif