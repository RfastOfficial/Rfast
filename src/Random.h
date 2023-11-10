#ifndef RANDOM_GENERATORS_H
#define RANDOM_GENERATORS_H

#include <limits>
#include <type_traits>
#include <chrono>
#include <vector>
#include <algorithm>

namespace Random {

	using std::numeric_limits;
	using real = std::false_type;
	using integer = std::true_type;
	using std::vector;
	using std::iota;

	namespace internal {

		static inline long long int get_current_nanoseconds(){
			return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		}

		class Integer_Core {
		public:
			// The next 3 lines are the requirements for UniformRandomBitGenerator.
			using result_type=uint32_t;
			static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
			static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }
		protected:
			struct pcg32_random_t {
				uint64_t state;
				uint64_t inc;
				pcg32_random_t(uint64_t init=get_current_nanoseconds()) : state(init), inc(init){}
			} rng;
			result_type pcg32_random_r()
			{
				uint64_t oldstate = rng.state;
		    // Advance internal state
				rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
		    // Calculate output function (XSH RR), uses old state for max ILP
				result_type xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
				result_type rot = oldstate >> 59u;
				return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
			}
		};
	}
	
	template<class T, bool replace = false>
	class uniform : public internal::Integer_Core{

		vector<size_t> indices;

		void remove_index(result_type i){
			this->indices[(size_t)i]=this->indices.back();
			this->indices.pop_back();
		}

	public:
		uniform(result_type max_bound=0){
			this->indices.resize(max_bound);
			iota(this->indices.begin(),this->indices.end(),0);
		}

		uniform(int32_t min_bound,int32_t max_bound){
			this->indices.resize(std::abs(max_bound-min_bound+1));
			iota(this->indices.begin(),this->indices.end(),min_bound);
		}

		result_type operator()(){
			auto index = (size_t)this->pcg32_random_r()%this->indices.size();
			auto res = this->indices[index];
			this->remove_index(index);
			return res;
		}

	};

	template<class T>
	class uniform<T,true> : public internal::Integer_Core{
	// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
	// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

		result_type min_bound,max_bound;

	public:
		uniform(result_type max_bound=1) : min_bound(1),max_bound(max_bound){}
		uniform(result_type min_bound,result_type max_bound) : min_bound(min_bound),max_bound(max_bound){}

		result_type operator()(){
			return this->pcg32_random_r()%this->max_bound+this->min_bound;
		}

	};

	template<>
	class uniform<real,false> : public internal::Integer_Core {

		const double min,max;

	public:
		uniform(const double min=0.0,const double max=1.0) : min(min), max(max){}

		inline double operator()(){
			return min + (this->pcg32_random_r() *  (max-min) / internal::Integer_Core::max() ) ;
		}

	};
}

#endif