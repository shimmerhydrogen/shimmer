#ifndef SHIMMER_TEQP_UTILITIES
#define SHIMMER_TEQP_UTILITIES

#include <list>
#include <ostream>
#include <vector>

namespace shimmer_teqp_utilities
{
  // *********************************************************
  template<typename collection_type>
  std::ostream& operator<<(std::ostream& out,
                           const collection_type& elements)
  {
    out<< "{";
    unsigned int i = 0;
    for (const auto& element : elements)
    {
      out<< (i != 0 ? ", " : "")<< element;
      i++;
    }
    out<< "}";

    return out;
  }
  // *********************************************************
  template<typename collection_type, typename value_type>
  std::vector<unsigned int> filter_components(const collection_type& mol_fracs,
                                              const value_type& tolerance)
  {
    std::list<unsigned int> filter;

    unsigned int i = 0;
    for (const auto& mol_frac : mol_fracs)
    {
      if (mol_frac > tolerance)
        filter.push_back(i);

      i++;
    }

    return std::vector<unsigned int>(filter.begin(), filter.end());
  }
  // *********************************************************
  template<typename collection_type, typename filter_type>
  collection_type filter_collection(const filter_type& filter,
                                    const collection_type& collection)
  {
    const auto filter_size = filter.size();
    using size_type = typename filter_type::size_type;

    size_type i = 0;
    collection_type filtered_collection(filter_size);
    for (const auto& f : filter)
    {
      filtered_collection[i++] = collection[f];
    }

    return filtered_collection;
  }
  // *********************************************************
}

#endif
