#ifndef __SHIMMER_GERG_UTILITIES_HPP__
#define __SHIMMER_GERG_UTILITIES_HPP__

#include <list>
#include <ostream>
#include <span>
#include <vector>

#include "shimmer_gerg_concepts.hpp"
#include "shimmer_gerg_data.hpp"

namespace shimmer_gerg
{
  namespace utilities
  {
    extern std::vector<std::string>  component_names;

    // *********************************************************
    template <concepts::IterableCollection collection>
    std::ostream& operator<<(std::ostream& out,
                             const collection& elements)
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
    // *************************************************************************
    template <concepts::IterableCollection collection1,
              concepts::IterableCollection collection2>
    bool operator==(const collection1& elements1,
                    const collection2& elements2)
    {
      assert(elements1.size() == elements2.size());
      auto view1 = std::span{elements1};
      auto view2 = std::span{elements2};

      auto it2 = elements2.begin();
      for (const auto& element1 : elements1)
      {
        if (element1 != (*it2))
          return false;
        it2++;
      }

      return true;
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
}

#endif
