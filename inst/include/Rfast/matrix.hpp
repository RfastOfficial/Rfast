
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"
#include "types.hpp"
#include "helpers.hpp"
#include "vector.hpp"
#include "FactorVector.hpp"

namespace Rfast
{
	using namespace Rcpp;
	using namespace arma;
	using namespace std;
	using namespace chrono;

	inline NumericMatrix transpose(NumericMatrix x, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.ncol(), n = x.nrow();
		NumericMatrix f = p == n ? clone(x) : NumericMatrix(p, n);
		if (p == n)
		{
			for (int i = 1; i < p; ++i)
			{
				for (int u = 0; u < i; ++u)
				{
					swap(f(u, i), f(i, u));
				}
			}
		}
		else
		{
			mat ff(f.begin(), p, n, false), xx(x.begin(), n, p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (int i = 0; i < p; ++i)
			{
				ff.row(i) = xx.col(i).t();
			}
		}
		return f;
	}

	inline mat transpose(mat x, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.n_cols, n = x.n_rows;
		mat f;
		if (p == n)
		{
			f = x;
			for (int i = 1; i < p; ++i)
			{
				for (int u = 0; u < i; ++u)
				{
					swap(f(u, i), f(i, u));
				}
			}
		}
		else
		{
			f = mat(p, n);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (int i = 0; i < p; ++i)
			{
				f.row(i) = x.col(i).t();
			}
		}
		return f;
	}

	inline NumericMatrix matrix_multiplication(NumericMatrix X, NumericMatrix Y, const bool tx = false, const bool ty = false, const unsigned int cores = get_num_of_threads())
	{
		int p, n;

		if (!tx)
		{
			n = X.nrow();
			p = ty ? Y.nrow() : Y.ncol();
		}
		else if (tx)
		{
			n = X.ncol();
			p = Y.ncol();
		}
		NumericMatrix C(n, p);
		mat CC(C.begin(), n, p, false), x(X.begin(), X.nrow(), X.ncol(), false), y(Y.begin(), Y.nrow(), Y.ncol(), false);
		colvec yi;

		if (!tx and !ty)
		{ // matrix multiplication
			mat xx = transpose(x);
			for (int i = 0; i < p; ++i)
			{
				yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					CC(j, i) = dot(xx.col(j), yi);
				}
			}
		}
		else if (tx)
		{ // crossprod
			for (int i = 0; i < p; ++i)
			{
				yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					CC(j, i) = dot(x.col(j), yi);
				}
			}
		}
		else
		{ // tcrossprod
			mat yy = transpose(y);
			mat xx = transpose(x);
			for (int i = 0; i < p; ++i)
			{
				yi = yy.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					CC(j, i) = dot(xx.col(j), yi);
				}
			}
		}
		return C;
	}

	inline mat matrix_multiplication(mat x, mat y, const bool tx = false, const bool ty = false, const unsigned int cores = get_num_of_threads())
	{
		int p = 0, n = 0;

		if (tx)
		{
			n = x.n_cols;
			p = x.n_cols;
		}
		else
		{
			n = x.n_rows;
			p = ty ? x.n_rows : x.n_cols;
		}
		mat C(n, p);
		colvec yi;

		if (!tx and !ty)
		{ // matrix multiplication
			mat xx = Rfast::transpose(x);
			for (int i = 0; i < p; ++i)
			{
				yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					C(j, i) = dot(xx.col(j), yi);
				}
			}
		}
		else if (tx)
		{ // crossprod
			for (int i = 0; i < p; ++i)
			{
				yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					C(j, i) = dot(x.col(j), yi);
				}
			}
		}
		else
		{ // tcrossprod
			mat yy = Rfast::transpose(y);
			mat xx = Rfast::transpose(x);
			for (int i = 0; i < p; ++i)
			{
				yi = yy.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int j = 0; j < n; ++j)
				{
					C(j, i) = dot(xx.col(j), yi);
				}
			}
		}
		return C;
	}

	inline NumericMatrix colSort(DataFrame x, const bool descend = false, const bool stable = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericMatrix F(x.nrows(), x.size());
		mat f(F.begin(), F.nrow(), F.ncol(), false);
		if (descend)
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						int i = s - x.begin();
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false, mgreater<bool, double, double>);
							break;
						case Type::Types::INT:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false, mgreater<bool, double, double>);
							break;
						case Type::Types::CHAR:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false, mgreater<bool, double, double>);
							break;
						case Type::Types::FACTOR:
							f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
				else
				{
					int i = 0;
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResult<colvec, std::stable_sort>(f, i, false,s, mgreater<bool, double, double>);
							break;
						case Type::Types::INT:
							setResult<colvec, std::stable_sort>(f, i, false,s, mgreater<bool, double, double>);
							break;
						case Type::Types::CHAR:
							setResult<colvec, std::stable_sort>(f, i, false,s, mgreater<bool, double, double>);
							break;
						case Type::Types::FACTOR:
							f.col(i++) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						int i = s - x.begin();
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i, false,mgreater<bool, double, double>);
							break;
						case Type::Types::INT:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i, false,mgreater<bool, double, double>);
							break;
						case Type::Types::CHAR:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i, false,mgreater<bool, double, double>);
							break;
						case Type::Types::FACTOR:
							f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
				else
				{
					int i = 0;
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResult<colvec, std::sort>(f, i, false, s, mgreater<bool, double, double>);
							break;
						case Type::Types::INT:
							setResult<colvec, std::sort>(f, i, false, s, mgreater<bool, double, double>);
							break;
						case Type::Types::CHAR:
							setResult<colvec, std::sort>(f, i, false, s, mgreater<bool, double, double>);
							break;
						case Type::Types::FACTOR:
							f.col(i++) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
			}
		}
		else
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						int i = s - x.begin();
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false);
							break;
						case Type::Types::INT:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false);
							break;
						case Type::Types::CHAR:
							setResultParallelSection<colvec, NumericVector, std::stable_sort>(f, s, i,false);
							break;
						case Type::Types::FACTOR:
							f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
				else
				{
					int i = 0;
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResult<colvec, std::stable_sort>(f, i++, false, s);
							break;
						case Type::Types::INT:
							setResult<colvec, std::stable_sort>(f, i++, false, s);
							break;
						case Type::Types::CHAR:
							setResult<colvec, std::stable_sort>(f, i++, false, s);
							break;
						case Type::Types::FACTOR:
							f.col(i++) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						int i = s - x.begin();
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i,false);
							break;
						case Type::Types::INT:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i,false);
							break;
						case Type::Types::CHAR:
							setResultParallelSection<colvec, NumericVector, std::sort>(f, s, i,false);
							break;
						case Type::Types::FACTOR:
							f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
				else
				{
					int i = 0;
					for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
					{
						switch (Type::type(s->get()))
						{
						case Type::Types::REAL:
							setResult<colvec, std::sort>(f, i++,false, s);
							break;
						case Type::Types::INT:
							setResult<colvec, std::sort>(f, i++,false, s);
							break;
						case Type::Types::CHAR:
							setResult<colvec, std::sort>(f, i++,false, s);
							break;
						case Type::Types::FACTOR:
							f.col(i++) = FactorVector::sort<colvec>(s->get(), descend);
							break;
						default:
							break;
						}
					}
				}
			}
		}
		colnames(F) = static_cast<CharacterVector>(x.names());
		return F;
	}

	inline NumericMatrix colSort(NumericMatrix X, const bool descend = false, const bool stable = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int n = X.nrow(), p = X.ncol();
		NumericMatrix F(n, p);
		mat f(F.begin(), n, p, false), x(X.begin(), n, p, false);
		if (descend)
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						stable_sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						stable_sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
			}
		}
		else
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						stable_sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						stable_sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
			}
		}
		return F;
	}

	inline mat colSort(mat x, const bool descend = false, const bool stable = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int n = x.n_rows, p = x.n_cols;
		mat f(n, p);
		if (descend)
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						stable_sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						stable_sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						sort(coli.begin(), coli.end(), greater<double>());
						f.col(i) = coli;
					}
				}
			}
		}
		else
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						stable_sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						stable_sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec coli = x.col(i);
						sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
				else
				{
					colvec coli(n);
					for (int i = 0; i < p; ++i)
					{
						coli = x.col(i);
						sort(coli.begin(), coli.end());
						f.col(i) = coli;
					}
				}
			}
		}
		return f;
	}

	inline NumericMatrix rowSort(NumericMatrix x, const bool descend = false, const bool stable = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int n = x.nrow(), p = x.ncol();
		NumericMatrix f(n, p);
		mat xx(x.begin(), n, p, false), ff(f.begin(), n, p, false);
		if (descend)
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = xx.row(i);
						stable_sort(rowi.begin(), rowi.end(), greater<double>());
						ff.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = xx.row(i);
						stable_sort(rowi.begin(), rowi.end(), greater<double>());
						ff.row(i) = rowi;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = xx.row(i);
						sort(rowi.begin(), rowi.end(), greater<double>());
						ff.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = xx.row(i);
						sort(rowi.begin(), rowi.end(), greater<double>());
						ff.row(i) = rowi;
					}
				}
			}
		}
		else
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = xx.row(i);
						stable_sort(rowi.begin(), rowi.end());
						ff.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = xx.row(i);
						stable_sort(rowi.begin(), rowi.end());
						ff.row(i) = rowi;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = xx.row(i);
						sort(rowi.begin(), rowi.end());
						ff.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = xx.row(i);
						sort(rowi.begin(), rowi.end());
						ff.row(i) = rowi;
					}
				}
			}
		}
		return f;
	}

	inline mat rowSort(mat x, const bool descend = false, const bool stable = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int n = x.n_rows, p = x.n_cols;
		mat f(n, p);
		if (descend)
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = x.row(i);
						stable_sort(rowi.begin(), rowi.end(), greater<double>());
						f.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = x.row(i);
						stable_sort(rowi.begin(), rowi.end(), greater<double>());
						f.row(i) = rowi;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = x.row(i);
						sort(rowi.begin(), rowi.end(), greater<double>());
						f.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = x.row(i);
						sort(rowi.begin(), rowi.end(), greater<double>());
						f.row(i) = rowi;
					}
				}
			}
		}
		else
		{
			if (stable)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = x.row(i);
						stable_sort(rowi.begin(), rowi.end());
						f.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = x.row(i);
						stable_sort(rowi.begin(), rowi.end());
						f.row(i) = rowi;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < n; ++i)
					{
						rowvec rowi = x.row(i);
						sort(rowi.begin(), rowi.end());
						f.row(i) = rowi;
					}
				}
				else
				{
					rowvec rowi(n);
					for (int i = 0; i < n; ++i)
					{
						rowi = x.row(i);
						sort(rowi.begin(), rowi.end());
						f.row(i) = rowi;
					}
				}
			}
		}
		return f;
	}

	inline bool is_symmetric(NumericMatrix x)
	{
		int ncl = x.ncol(), i, j;
		for (i = 1; i < ncl; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				if (x(j, i) != x(i, j))
				{
					return false;
				}
			}
		}
		return true;
	}

	inline bool is_symmetric(mat x)
	{
		int ncl = x.n_cols, i, j;
		for (i = 1; i < ncl; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				if (x(j, i) != x(i, j))
				{
					return false;
				}
			}
		}
		return true;
	}

	inline NumericVector colMedian(DataFrame &x, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector f(x.size());
		colvec ff(f.begin(),f.size(),false);
		if (parallel)
		{
			colvec ff(f.begin(), f.size(), false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				int i = s - x.begin();
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					setResultParallelSection<colvec, NumericVector, med_helper<colvec>>(ff, s, i, na_rm);
					break;
				case Type::Types::INT:
					setResultParallelSection<colvec, NumericVector, med_helper<colvec>>(ff, s, i, na_rm);
					break;
				case Type::Types::CHAR:
					setResultParallelSection<colvec, NumericVector, med_helper<colvec>>(ff, s, i, na_rm);
					break;
				case Type::Types::FACTOR:
					//f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		else
		{
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				int i = s - x.begin();
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					setResult<NumericVector, med_helper<colvec>>(ff, i, na_rm, s);
					break;
				case Type::Types::INT:
					setResult<NumericVector, med_helper<colvec>>(ff, i, na_rm, s);
					break;
				case Type::Types::CHAR:
					setResult<NumericVector, med_helper<colvec>>(ff, i, na_rm, s);
					break;
				case Type::Types::FACTOR:
					//f[i] = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		f.names() = static_cast<CharacterVector>(x.names());
		return f;
	}

	inline NumericVector colMedian(NumericMatrix &x, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.ncol();
		NumericVector F(p);
		if (na_rm)
		{
			if (parallel)
			{
				mat xx(x.begin(), x.nrow(), p, false);
				colvec ff(F.begin(), p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int i = 0; i < p; ++i)
				{
					colvec tmp = xx.col(i);
					ff[i] = med_helper<colvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
			else
			{
				NumericVector tmp(x.nrow());
				for (int i = 0; i < p; ++i)
				{
					tmp = x.column(i);
					F[i] = med_helper<NumericVector>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
		}
		else
		{
			const int step = x.nrow(), middle = step / 2 - 1;
			if (step % 2 == 0)
			{
				if (parallel)
				{
					mat xx(x.begin(), step, p, false);
					colvec ff(F.begin(), p, false), tmpp(step);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec tmpp = xx.col(i);
						nth_element(tmpp.begin(), tmpp.begin() + middle, tmpp.end());
						ff[i] = (tmpp[middle] + *(min_element(tmpp.begin() + middle + 1, tmpp.end()))) / 2.0;
					}
				}
				else
				{
					NumericVector tmp(step);
					for (int i = 0; i < p; ++i)
					{
						tmp = x.column(i);
						nth_element(tmp.begin(), tmp.begin() + middle, tmp.end());
						F[i] = (tmp[middle] + *(min_element(tmp.begin() + middle + 1, tmp.end()))) / 2.0;
					}
				}
			}
			else
			{
				if (parallel)
				{
					mat xx(x.begin(), step, p, false);
					colvec ff(F.begin(), p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec tmpp = xx.col(i);
						nth_element(tmpp.begin(), tmpp.begin() + middle + 1, tmpp.end());
						ff[i] = tmpp[middle + 1];
					}
				}
				else
				{
					NumericVector tmp(step);
					for (int i = 0; i < p; ++i)
					{
						tmp = x.column(i);
						nth_element(tmp.begin(), tmp.begin() + middle + 1, tmp.end());
						F[i] = tmp[middle + 1];
					}
				}
			}
		}
		return F;
	}

	inline rowvec colMedian(mat &x, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.n_cols;
		rowvec F(p);
		if (na_rm)
		{
			if (parallel)
			{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int i = 0; i < p; ++i)
				{
					colvec tmp = x.col(i);
					F[i] = med_helper<colvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
			else
			{
				colvec tmp(x.n_rows);
				for (int i = 0; i < p; ++i)
				{
					tmp = x.col(i);
					F[i] = med_helper<colvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
		}
		else
		{
			const int step = x.n_rows, middle = step / 2 - 1;
			if (step % 2 == 0)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec tmp = x.col(i);
						nth_element(tmp.begin(), tmp.begin() + middle, tmp.end());
						F[i] = (tmp[middle] + *(min_element(tmp.begin() + middle + 1, tmp.end()))) / 2.0;
					}
				}
				else
				{
					colvec tmp(step);
					for (int i = 0; i < p; ++i)
					{
						tmp = x.col(i);
						nth_element(tmp.begin(), tmp.begin() + middle, tmp.end());
						F[i] = (tmp[middle] + *(min_element(tmp.begin() + middle + 1, tmp.end()))) / 2.0;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						colvec tmp = x.col(i);
						nth_element(tmp.begin(), tmp.begin() + middle + 1, tmp.end());
						F[i] = tmp[middle + 1];
					}
				}
				else
				{
					colvec tmp(step);
					for (int i = 0; i < p; ++i)
					{
						tmp = x.col(i);
						nth_element(tmp.begin(), tmp.begin() + middle + 1, tmp.end());
						F[i] = tmp[middle + 1];
					}
				}
			}
		}
		return F;
	}

	inline NumericVector rowMedian(NumericMatrix x, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.nrow();
		NumericVector F(p);
		if (na_rm)
		{
			if (parallel)
			{
				mat xx(x.begin(), x.nrow(), p, false);
				colvec ff(F.begin(), p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int i = 0; i < p; ++i)
				{
					rowvec tmp = xx.row(i);
					ff[i] = med_helper<rowvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
			else
			{
				NumericVector tmp(x.ncol());
				for (int i = 0; i < p; ++i)
				{
					tmp = x.row(i);
					F[i] = med_helper<rowvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
		}
		else
		{
			const int sz = x.ncol(), middle = sz / 2 - 1;
			if (sz % 2 == 0)
			{
				if (parallel)
				{
					mat xx(x.begin(), p, sz, false);
					colvec ff(F.begin(), p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						rowvec rowii = xx.row(i);
						nth_element(rowii.begin(), rowii.begin() + middle, rowii.end());
						ff[i] = (rowii[middle] + *(min_element(rowii.begin() + middle + 1, rowii.end()))) / 2.0;
					}
				}
				else
				{
					NumericVector rowi(sz);
					for (int i = 0; i < p; ++i)
					{
						rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = (rowi[middle] + *(min_element(rowi.begin() + middle + 1, rowi.end()))) / 2.0;
					}
				}
			}
			else
			{
				if (parallel)
				{
					mat xx(x.begin(), p, sz, false);
					colvec ff(F.begin(), p, false);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						rowvec rowii = xx.row(i);
						nth_element(rowii.begin(), rowii.begin() + middle, rowii.end());
						ff[i] = rowii[middle + 1];
					}
				}
				else
				{
					NumericVector rowi(sz);
					for (int i = 0; i < p; ++i)
					{
						rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = rowi[middle + 1];
					}
				}
			}
		}
		return F;
	}

	inline colvec rowMedian(mat &x, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		const int p = x.n_rows;
		colvec F(p);
		if (na_rm)
		{
			if (parallel)
			{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
				for (int i = 0; i < p; ++i)
				{
					rowvec tmp = x.row(i);
					F[i] = med_helper<rowvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
			else
			{
				rowvec tmp(x.n_cols);
				for (int i = 0; i < p; ++i)
				{
					tmp = x.row(i);
					F[i] = med_helper<rowvec>(tmp.begin(), tmp.begin() + (int)(std::remove_if(tmp.begin(), tmp.end(), R_IsNA) - tmp.begin()));
				}
			}
		}
		else
		{
			const int sz = x.n_cols, middle = sz / 2 - 1;
			if (sz % 2 == 0)
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						rowvec rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = (rowi[middle] + *(min_element(rowi.begin() + middle + 1, rowi.end()))) / 2.0;
					}
				}
				else
				{
					rowvec rowi(sz);
					for (int i = 0; i < p; ++i)
					{
						rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = (rowi[middle] + *(min_element(rowi.begin() + middle + 1, rowi.end()))) / 2.0;
					}
				}
			}
			else
			{
				if (parallel)
				{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
					for (int i = 0; i < p; ++i)
					{
						rowvec rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = rowi[middle + 1];
					}
				}
				else
				{
					rowvec rowi(sz);
					for (int i = 0; i < p; ++i)
					{
						rowi = x.row(i);
						nth_element(rowi.begin(), rowi.begin() + middle, rowi.end());
						F[i] = rowi[middle + 1];
					}
				}
			}
		}
		return F;
	}

	inline NumericVector colVars(NumericMatrix x, const bool std = false, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat xx(x.begin(), x.nrow(), x.ncol(), false);
		NumericVector f(xx.n_cols);
		rowvec ff(f.begin(), f.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < xx.n_cols; ++i)
			{
				ff[i] = Rfast::var<colvec>(xx.col(i), std, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < xx.n_cols; ++i)
			{
				ff[i] = Rfast::var<colvec>(xx.col(i), std, na_rm);
			}
		}
		return f;
	}

	inline NumericVector colVars(DataFrame x, const bool std = false, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector f(x.size());
		colvec ff(f.begin(),f.size(),false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::INT:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::CHAR:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::FACTOR:
					//f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		else
		{
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::INT:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::CHAR:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::var<colvec>, std, na_rm);
					break;
				case Type::Types::FACTOR:
					//f[i] = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		f.names() = static_cast<CharacterVector>(x.names());
		return f;
	}

	inline NumericVector rowVars(NumericMatrix x, const bool std = false, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat xx(x.begin(), x.nrow(), x.ncol(), false);
		NumericVector f(xx.n_rows);
		colvec ff(f.begin(), f.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < xx.n_rows; ++i)
			{
				ff[i] = Rfast::var<rowvec>(xx.row(i), std, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < xx.n_rows; ++i)
			{
				ff[i] = Rfast::var<rowvec>(xx.row(i), std, na_rm);
			}
		}
		return f;
	}

	inline rowvec colVars(mat x, const bool std = false, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		rowvec f(x.n_cols);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < x.n_cols; ++i)
			{
				f[i] = Rfast::var<colvec>(x.col(i), std, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < x.n_cols; ++i)
			{
				f[i] = Rfast::var<colvec>(x.col(i), std, na_rm);
			}
		}
		return f;
	}

	inline colvec rowVars(mat x, const bool std = false, const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		colvec f(x.n_rows);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < x.n_rows; ++i)
			{
				f[i] = Rfast::var<rowvec>(x.row(i), std, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < x.n_rows; ++i)
			{
				f[i] = Rfast::var<rowvec>(x.row(i), std, na_rm);
			}
		}
		return f;
	}

	inline NumericVector colMads(DataFrame x, const string method = "median", const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector f(x.size());
		colvec ff(f.begin(),f.size(),false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::INT:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::CHAR:
					ff[s-x.begin()] = setResultParallelSection<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::FACTOR:
					//f.col(i) = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		else
		{
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::INT:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::CHAR:
					ff[s-x.begin()] = singleIteratorWithoutCopy<colvec, NumericVector>(s, Rfast::mad<colvec>, method, na_rm);
					break;
				case Type::Types::FACTOR:
					//f[i] = FactorVector::sort<colvec>(s->get(), descend);
					break;
				default:
					break;
				}
			}
		}
		f.names() = static_cast<CharacterVector>(x.names());
		return f;
	}

	inline NumericVector colMads(NumericMatrix x, const string method = "median", const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat xx(x.begin(), x.nrow(), x.ncol(), false);
		NumericVector f(xx.n_cols);
		rowvec ff(f.begin(), f.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < xx.n_cols; ++i)
			{
				ff[i] = Rfast::mad<colvec>(xx.col(i), method, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < xx.n_cols; ++i)
			{
				ff[i] = Rfast::mad<colvec>(xx.col(i), method, na_rm);
			}
		}
		return f;
	}

	inline NumericVector rowMads(NumericMatrix x, const string method = "median", const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat xx(x.begin(), x.nrow(), x.ncol(), false);
		NumericVector f(xx.n_rows);
		colvec ff(f.begin(), f.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < xx.n_rows; ++i)
			{
				ff[i] = Rfast::mad<rowvec>(xx.row(i), method, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < xx.n_rows; ++i)
			{
				ff[i] = Rfast::mad<rowvec>(xx.row(i), method, na_rm);
			}
		}
		return f;
	}

	inline rowvec colMads(mat x, const string method = "median", const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		rowvec f(x.n_cols);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < x.n_cols; ++i)
			{
				f[i] = Rfast::mad<colvec>(x.col(i), method, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < x.n_cols; ++i)
			{
				f[i] = Rfast::mad<colvec>(x.col(i), method, na_rm);
			}
		}
		return f;
	}

	inline colvec rowMads(mat x, const string method = "median", const bool na_rm = false, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		colvec f(x.n_rows);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (unsigned int i = 0; i < x.n_rows; ++i)
			{
				f[i] = Rfast::mad<rowvec>(x.row(i), method, na_rm);
			}
		}
		else
		{
			for (unsigned int i = 0; i < x.n_rows; ++i)
			{
				f[i] = Rfast::mad<rowvec>(x.row(i), method, na_rm);
			}
		}
		return f;
	}

	template <class Engine = std::default_random_engine>
	inline DataFrame colShuffle(DataFrame x, Engine engine = Engine())
	{
		const int n = x.size();
		seed_seq seq{get_current_nanoseconds()};
		std::vector<long long unsigned int> seeds(n);
		seq.generate(seeds.begin(), seeds.end());
		List f(n);
		// if(parallel){
		// 	colvec ff(f.begin(),f.size(),false);
		// 	#pragma omp parallel for num_threads(cores)
		// 	for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
		// 		colvec y;
		// 		int i;
		// 		#pragma omp critical
		// 		{
		// 			engine.seed(seeds[i]);
		// 			NumericVector yy;
		// 			yy=*s;
		// 			y = colvec(yy.begin(),yy.size(),false);
		// 			i = s-x.begin();
		// 		}
		// 		ff[i]=Rfast::shuffle<NumericVector>(y,engine);
		// 	}
		// }else{
		int i = 0;
		for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
		{
			engine.seed(seeds[i]);
			switch (Type::type(s->get()))
			{
			case Type::Types::REAL:
				setResult<NumericVector>(f, i, s, Rfast::shuffle<colvec>, engine);
				break;
			case Type::Types::INT:
				setResult<NumericVector>(f, i, s, Rfast::shuffle<colvec>, engine);
				break;
			case Type::Types::CHAR:
				setResult<NumericVector>(f, i, s, Rfast::shuffle<colvec>, engine);
				break;
			case Type::Types::FACTOR:
				//f[i] = FactorVector::sort<colvec>(s->get(), descend);
				break;
			default:
				break;
			}
		}
		f.names() = static_cast<CharacterVector>(x.names());
		return f;
	}

	template <class Engine = std::default_random_engine>
	inline NumericMatrix colShuffle(NumericMatrix x, Engine engine = Engine())
	{
		const int n = x.ncol();
		seed_seq seq{get_current_nanoseconds()};
		std::vector<long long unsigned int> seeds(n);
		seq.generate(seeds.begin(), seeds.end());
		NumericMatrix y(x.nrow(), n);
		for (int i = 0; i < n; ++i)
		{
			engine.seed(seeds[i]);
			y.column(i) = Rfast::shuffle<NumericVector>(x.column(i), engine);
		}
		return y;
	}

	template <class Engine = std::default_random_engine>
	NumericMatrix rowShuffle(NumericMatrix x, Engine engine = Engine())
	{
		const int n = x.ncol();
		seed_seq seq{get_current_nanoseconds()};
		std::vector<long long unsigned int> seeds(n);
		seq.generate(seeds.begin(), seeds.end());
		NumericMatrix y(x.nrow(), n);
		for (int i = 0; i < n; ++i)
		{
			engine.seed(seeds[i]);
			y.row(i) = Rfast::shuffle<NumericVector>(x.row(i), engine);
		}
		return y;
	}

	template <class Engine = std::default_random_engine>
	mat colShuffle(mat x, Engine engine = Engine())
	{
		const int n = x.n_cols;
		seed_seq seq{get_current_nanoseconds()};
		std::vector<long long unsigned int> seeds(n);
		seq.generate(seeds.begin(), seeds.end());
		mat y(x.n_rows, n);
		for (int i = 0; i < n; ++i)
		{
			engine.seed(seeds[i]);
			y.col(i) = Rfast::shuffle<colvec>(x.col(i), engine);
		}
		return y;
	}

	template <class Engine = std::default_random_engine>
	mat rowShuffle(mat x, Engine engine = Engine())
	{
		const int n = x.n_rows;
		seed_seq seq{get_current_nanoseconds()};
		std::vector<long long unsigned int> seeds(n);
		seq.generate(seeds.begin(), seeds.end());
		mat y(n, x.n_cols);
		for (int i = 0; i < n; ++i)
		{
			engine.seed(seeds[i]);
			y.row(i) = Rfast::shuffle<rowvec>(x.row(i), engine);
		}
		return y;
	}

	inline NumericVector colMaxs(DataFrame x, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector F(x.size());
		colvec f(F.begin(), F.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<colvec, NumericVector, std::max_element>(s);
					break;
				case Type::Types::INT:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<icolvec, IntegerVector, std::max_element>(s);
					break;
				case Type::Types::CHAR:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<icolvec, IntegerVector, std::max_element>(s);
					break;
				case Type::Types::FACTOR:
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					f[s - x.begin()] = FactorVector(s->get()).maxIndex();
				}
				break;
				default:
					break;
				}
			}
		}
		else
		{
			int i = 0;
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f[i++] = singleIteratorWithoutCopy<colvec, NumericVector, std::max_element>(s);
					break;
				case Type::Types::INT:
					f[i++] = singleIteratorWithoutCopy<icolvec, IntegerVector, std::max_element>(s);
					break;
				case Type::Types::CHAR:
					f[i++] = singleIteratorWithoutCopy<icolvec, IntegerVector, std::max_element>(s);
					break;
				case Type::Types::FACTOR:
					f[i++] = FactorVector(s->get()).maxIndex();
					break;
				default:
					break;
				}
			}
		}
		colnames(F) = static_cast<CharacterVector>(x.names());
		return F;
	}

	inline NumericVector colMins(DataFrame x, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector F(x.size());
		colvec f(F.begin(), F.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<colvec, NumericVector, std::min_element>(s);
					break;
				case Type::Types::INT:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<icolvec, IntegerVector, std::min_element>(s);
					break;
				case Type::Types::CHAR:
					f[s - x.begin()] = parallelSingleIteratorWithoutCopy<icolvec, IntegerVector, std::min_element>(s);
					break;
				case Type::Types::FACTOR:
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					f[s - x.begin()] = FactorVector(s->get()).maxIndex();
				}
				break;
				default:
					break;
				}
			}
		}
		else
		{
			int i = 0;
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f[i++] = singleIteratorWithoutCopy<colvec, NumericVector, std::min_element>(s);
					break;
				case Type::Types::INT:
					f[i++] = singleIteratorWithoutCopy<icolvec, IntegerVector, std::min_element>(s);
					break;
				case Type::Types::CHAR:
					f[i++] = singleIteratorWithoutCopy<icolvec, IntegerVector, std::min_element>(s);
					break;
				case Type::Types::FACTOR:
					f[i++] = FactorVector(s->get()).maxIndex();
					break;
				default:
					break;
				}
			}
		}
		colnames(F) = static_cast<CharacterVector>(x.names());
		return F;
	}

	inline NumericMatrix colMinsMaxs(DataFrame x, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericMatrix F(2, x.size());
		mat f(F.begin(), 2, F.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f.col(s - x.begin()) = parallelSingleIteratorWithoutCopy<colvec, colvec, NumericVector, std::minmax_element>(s);
					break;
				case Type::Types::INT:
					f.col(s - x.begin()) = parallelSingleIteratorWithoutCopy<colvec, icolvec, IntegerVector, std::minmax_element>(s);
					break;
				case Type::Types::CHAR:
					f.col(s - x.begin()) = parallelSingleIteratorWithoutCopy<colvec, icolvec, IntegerVector, std::minmax_element>(s);
					break;
				case Type::Types::FACTOR:
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					f.col(s - x.begin()) = FactorVector(s->get()).minmaxIndex<colvec>();
				}
				break;
				default:
					break;
				}
			}
		}
		else
		{
			int i = 0;
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				switch (Type::type(s->get()))
				{
				case Type::Types::REAL:
					f.col(i++) = singleIteratorWithoutCopy<colvec, colvec, NumericVector, std::minmax_element>(s);
					break;
				case Type::Types::INT:
					f.col(i++) = singleIteratorWithoutCopy<colvec, icolvec, IntegerVector, std::minmax_element>(s);
					break;
				case Type::Types::CHAR:
					f.col(i++) = singleIteratorWithoutCopy<colvec, icolvec, IntegerVector, std::minmax_element>(s);
					break;
				case Type::Types::FACTOR:
					f.col(i++) = FactorVector(s->get()).minmaxIndex<colvec>();
					break;
				default:
					break;
				}
			}
		}
		colnames(F) = static_cast<CharacterVector>(x.names());
		rownames(F) = CharacterVector::create("min", "max");
		return F;
	}
}

#endif
