
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

namespace Rfast
{

	class FactorVector : public IntegerVector
	{
	private:
		CharacterVector levels;

		friend iterator;

	public:
		FactorVector(const FactorVector &f) : IntegerVector(wrap(f))
		{
			this->levels = f.levels;
		}

		template <typename Proxy>
		FactorVector(const GenericProxy<Proxy> &proxy) : IntegerVector(proxy), levels(attr("levels")) {}

		FactorVector(SEXP x) : IntegerVector(x), levels(attr("levels")) {}

		template <class T>
		FactorVector(T x) : IntegerVector(x), levels(attr("levels")) {}

		FactorVector() {}

		/*void operator=(const FactorVector rhs)
		 {
		 Rcout<<__LINE__<<endl;
		 this->levels = rhs.levels;
		 }*/

		void setLevels(CharacterVector levels)
		{
			this->levels = levels;
		}

		CharacterVector getLevels() const
		{
			return this->levels;
		}

		class iterator
		{
			const FactorVector &base;
			size_t index;

		public:
			iterator(const FactorVector &base) : base(base), index(0) {}

			CharacterVector::const_Proxy operator*()
			{
				return base.levels[base[index]];
			}

			IntegerVector::stored_type level() const
			{
				return base[index];
			}

			iterator &operator++()
			{
				++this->index;
				return *this;
			}

			iterator operator++(int)
			{
				iterator res(base);
				res.index = this->index;
				return res;
			}

			iterator &operator+(size_t n)
			{
				this->index += n;
				return *this;
			}

			bool operator==(iterator i)
			{
				return this->operator*() == *i;
			}

			bool operator!=(iterator i)
			{
				return this->operator*() != *i;
			}

			iterator &operator=(iterator &i)
			{
				this->index = i.index;
				return *this;
			}

			inline IntegerVector::iterator begin_ptr() const
			{
				return static_cast<IntegerVector>(base).begin();
			}

			const FactorVector &parent() const
			{
				return this->base;
			}
		};

		inline iterator begin() const
		{
			return iterator(*this);
		}

		inline iterator end() const
		{
			return iterator(*this) + size();
		}

		CharacterVector::Proxy max()
		{
			return levels[levels.size() - 1];
		}

		CharacterVector::Proxy min()
		{
			return levels[0];
		}

		size_t maxIndex()
		{
			return levels.size();
		}

		size_t minIndex()
		{
			return 1;
		}

		template <class T>
		T minmaxIndex()
		{
			return {static_cast<typename T::value_type>(minIndex()), static_cast<typename T::value_type>(maxIndex())};
		}

		template <class T>
		T minmax()
		{
			return {static_cast<typename T::value_type>(min()), static_cast<typename T::value_type>(max())};
		}

		template <class T>
		T sort(const bool descend = false)
		{
			icolvec x, inds;
			IntegerVector I(*this);
			x = icolvec(I.begin(), I.size());
			inds = Tabulate<icolvec, icolvec>(x, x.n_elem);
			T res(x.size());
			int start = 0;

			if (descend)
			{
				start = res.size() - 1;
				for (size_t i = 0; i < inds.size(); ++i)
				{
					int v = inds[i];
					if (v > 0)
					{
						const int end = start - v;
						for (; start > end; --start)
						{
							res[start] = i + 1;
						}
					}
					else
					{
						break;
					}
				}
			}
			else
			{
				start = 0;
				for (size_t i = 0; i < inds.size(); ++i)
				{
					int v = inds[i];
					if (v > 0)
					{
						const int end = start + v;
						for (; start < end; ++start)
						{
							res[start] = i + 1;
						}
					}
					else
					{
						break;
					}
				}
			}
			return res;
		}

		template <class T>
		static T sort(SEXP xx, const bool descend = false)
		{
			icolvec x, inds;
#pragma omp critical
			{
				IntegerVector I(xx);
				x = icolvec(I.begin(), I.size(), false);
				inds = Tabulate<icolvec, icolvec>(x, x.n_elem);
			}
			T res(x.size());
			int start = 0;

			if (descend)
			{
				start = res.size() - 1;
				for (size_t i = 0; i < inds.size(); ++i)
				{
					int v = inds[i];
					if (v > 0)
					{
						const int end = start - v;
						for (; start > end; --start)
						{
							res[start] = i + 1;
						}
					}
					else
					{
						break;
					}
				}
			}
			else
			{
				start = 0;
				for (size_t i = 0; i < inds.size(); ++i)
				{
					int v = inds[i];
					if (v > 0)
					{
						const int end = start + v;
						for (; start < end; ++start)
						{
							res[start] = i + 1;
						}
					}
					else
					{
						break;
					}
				}
			}
			return res;
		}
	};

}
#endif