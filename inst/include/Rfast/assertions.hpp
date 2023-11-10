

#ifndef ASSERTIONS_H
#define ASSERTIONS_H

namespace Assertion
{

#include <utility>
#include <string_view>

    using namespace std;

    template <typename, typename T = void>
    struct is_iterable : false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!is_same<T, void>::value, "The template class must be iterable.");
        }
    };

    template <typename T>
    struct is_iterable<T,
                       void_t<decltype(declval<T>().begin()), decltype(declval<T>().end())>> : true_type
    {
        constexpr static void check_concept() {}
    };

    template <typename, typename T = void>
    struct has_size : false_type
    {
        constexpr static void check_concept()
        {
            static_assert(!is_same<T, void>::value, "The template class must provide a function named size.");
        }
    };

    template <typename T>
    struct has_size<T, void_t<decltype(declval<T>().size())>> : true_type
    {
        constexpr static void check_concept() {}
    };
}

#endif