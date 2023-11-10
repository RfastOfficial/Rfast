#include <Rcpp.h>
using namespace Rcpp;
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