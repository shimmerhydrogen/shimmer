#ifndef __SHIMMER_TEQP_CONCEPTS_HPP__
#define __SHIMMER_TEQP_CONCEPTS_HPP__

#include <cassert>
#include <string_view>
#include <utility>
#include <span>

namespace shimmer_teqp
{
  namespace concepts
  {
    // *************************************************************************
    template<class T>
    concept IsStringLike = std::is_convertible_v<T, std::string_view>;
    // *************************************************************************
    template <typename T>
    concept IsNotStringLike = not IsStringLike<T>;
    // *************************************************************************
    template <typename collection>
    concept ResizableCollection = requires(collection t)
    { t.resize(0); };
    // *************************************************************************
    template <typename collection>
    concept FixedCollection = not ResizableCollection<collection>;
    // *************************************************************************
    template <typename collection>
    concept IsSpan = std::is_same<collection, std::span<typename collection::value_type>>::value ||
    std::is_same<collection, std::span<const typename collection::value_type>>::value;
    // *************************************************************************
    template <typename collection>
    concept IsNotSpan = not IsSpan<collection>;
    // *************************************************************************
    template <class collection>
    concept IterableCollection = IsNotStringLike<collection> && requires(collection a)
    { { a.begin() } -> std::same_as<typename collection::iterator>;
    { a.end() } -> std::same_as<typename collection::iterator>;
    { a.size() } -> std::same_as<typename collection::size_type>; };
    // *************************************************************************
    template <typename matrix>
    concept IsMatrix = requires(matrix m)
    { { m.cols() } -> std::same_as<typename matrix::size_type>;
    { m.rows() } -> std::same_as<typename matrix::size_type>; };
    // *************************************************************************
    template <typename expression>
    concept IsExpression = requires(expression e)
    { { expression::IsExpression() } -> std::same_as<bool>; };
    // *************************************************************************
  }
}
#endif // __CONCEPTS_HPP__
