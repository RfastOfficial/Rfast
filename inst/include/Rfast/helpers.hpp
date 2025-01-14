
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

template <class T, class HELPER, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function>
typename T::value_type parallelSingleIteratorWithoutCopy(DataFrame::iterator s)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = T(yy.begin(), yy.size(), false);
#ifdef _OPENMP
	}
#endif
	return *Function(y.begin(), y.end());
}

template <class T, class Helper, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function>
typename T::value_type singleIteratorWithoutCopy(DataFrame::iterator c)
{
	Helper h(*c);
	T y(h.begin(), h.size(), false);
	return *Function(y.begin(), y.end());
}

template <class T, class Helper, class Function, class ...Args>
typename T::value_type singleIteratorWithoutCopy(DataFrame::iterator c, Function func, Args... args)
{
	Helper h(*c);
	T y(h.begin(), h.size(), false);
	return func(y,args...);
}

template <class T, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function>
void singleIteratorWithoutCopy(List &f, const int i, DataFrame::iterator c)
{
	T y(*c);
	f[i] = T::create(*Function(y.begin(), y.end()));
}

template <class T, class Helper, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function>
void singleIteratorWithoutCopy(mat &f, const int i, DataFrame::iterator c)
{
	Helper h(*c);
	T y(h.begin(), h.size(), false);
	f[i] = *Function(y.begin(), y.end());
}

template <class T, Mfunction<void, typename T::iterator, typename T::iterator> Function>
void setResult(List &f, const int i, DataFrame::iterator c)
{
	T y = clone(as<T>(*c));
	Function(y.begin(), y.end());
	f[i] = T(y.begin(), y.end());
}

template <class T, class Function, class ...Args>
void setResult(List &f, const int i, DataFrame::iterator c, Function func, Args ...args)
{
	T y = clone(as<T>(*c));
	func(y, args...);
	f[i] = T(y.begin(), y.end());
}

template <class T, Mfunction<typename std::remove_reference<typename T::value_type>::type, typename T::iterator, typename T::iterator> Function>
void setResult(colvec &f, const int i,const bool na_rm, DataFrame::iterator c)
{
	T y = clone(as<T>(*c));
	f[i] = na_rm ? Function(y.begin(), y.end()) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()));
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
		  class F>
void setResult(mat &f, const int i,const bool na_rm, DataFrame::iterator c, F cmp)
{
	T y = as<T>(*c);
	na_rm ? Function(y.begin(), y.end(), cmp) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()), cmp);
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
		  class F>
void setResultParallelSection(mat &f, DataFrame::iterator s, const int i,const bool na_rm, F cmp)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = as<T>(yy);
#ifdef _OPENMP
	}
#endif
	na_rm ? Function(y.begin(), y.end(), cmp) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()), cmp);
	f.col(i) = y;
}

template <class T, class HELPER,Mfunction<typename T::value_type,typename T::iterator,typename T::iterator> Function>
void setResultParallelSection(colvec &f, DataFrame::iterator s, const int i,const bool na_rm)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = as<T>(yy);
#ifdef _OPENMP
	}
#endif
	f[i] = na_rm ? Function(y.begin(), y.end()) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()));
}

template <class T, class HELPER, class Function,class ...Args>
typename T::value_type setResultParallelSection(DataFrame::iterator s, Function func, Args... args)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = as<T>(yy);
#ifdef _OPENMP
	}
#endif
	return func(y,args...);
}

template <class T, class HELPER, Mfunction<void, typename T::iterator, typename T::iterator> Function>
void setResultParallelSection(mat &f, DataFrame::iterator s, const int i,const bool na_rm)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = as<T>(yy);
#ifdef _OPENMP
	}
#endif
	na_rm ? Function(y.begin(), y.end()) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()));
	f.col(i) = y;
}

template <class T, Mfunction<void, typename T::iterator, typename T::iterator> Function>
void setResult(mat &f, const int i,const bool na_rm, DataFrame::iterator c)
{
	T y = as<T>(*c);
	na_rm ? Function(y.begin(), y.end()) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()));
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
		  class F>
void setResult(List &f, const int i,const bool na_rm, DataFrame::iterator c, F cmp)
{
	T y = clone(as<T>(*c));
	na_rm ? Function(y.begin(), y.end(), cmp) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()), cmp);
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
		  class F>
void setResultParallelSection(List &f,const bool na_rm, DataFrame::iterator s, F cmp)
{
	T y;
	int i;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy = *s;
		y = as<T>(yy);
		i = s - f.begin();
#ifdef _OPENMP
	}
#endif
	na_rm ? Function(y.begin(), y.end(), cmp) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()), cmp);
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		f[i] = HELPER(y.begin(), y.end());
#ifdef _OPENMP
	}
#endif
}

template <class T, class HELPER, Mfunction<void, typename T::iterator, typename T::iterator> Function>
void setResultParallelSection(List &f,const bool na_rm, DataFrame::iterator s)
{
	T y;
	int i;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy = *s;
		y = as<T>(yy);
		i = s - f.begin();
#ifdef _OPENMP
	}
#endif
	na_rm ? Function(y.begin(), y.end()) : Function(y.begin(), y.begin()+(int)(std::remove_if(y.begin(), y.end(), R_IsNA) - y.begin()));

#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		f[i] = HELPER(y.begin(), y.end());
#ifdef _OPENMP
	}
#endif
}

template <class T, class HELPER, Mfunction<typename T::iterator, typename T::iterator, typename T::iterator> Function>
void parallelSingleIteratorWithoutCopy(List &f, DataFrame::iterator s)
{
	T y;
	int i;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = as<T>(yy);
		i = s - f.begin();
#ifdef _OPENMP
	}
#endif
	auto value = *Function(y.begin(), y.end());

#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		f[i] = HELPER::create(value);
#ifdef _OPENMP
	}
#endif
}

template <class RET, class T, class HELPER, Mfunction<std::pair<typename T::iterator, typename T::iterator>, typename T::iterator, typename T::iterator> Function>
RET parallelSingleIteratorWithoutCopy(DataFrame::iterator s)
{
	T y;
#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		HELPER yy(*s);
		y = T(yy.begin(), yy.size(), false);
#ifdef _OPENMP
	}
#endif
	auto v = Function(y.begin(), y.end());
	return {static_cast<typename RET::value_type>(*v.first), static_cast<typename RET::value_type>(*v.second)};
}

template <class RET, class T, class Helper, Mfunction<std::pair<typename T::iterator, typename T::iterator>, typename T::iterator, typename T::iterator> Function>
RET singleIteratorWithoutCopy(DataFrame::iterator c)
{
	Helper h(*c);
	T y(h.begin(), h.size(), false);
	auto v = Function(y.begin(), y.end());
	return {static_cast<typename RET::value_type>(*v.first), static_cast<typename RET::value_type>(*v.second)};
}

inline long long int get_current_nanoseconds()
{
	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}
#endif