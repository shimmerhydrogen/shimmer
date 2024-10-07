#ifndef __SHIMMER_TEQP_UTILITIES_HPP__
#define __SHIMMER_TEQP_UTILITIES_HPP__

#include <list>
#include <ostream>
#include <span>
#include <vector>

#include "shimmer_teqp_concepts.hpp"
#include "shimmer_gerg_data.hpp"

namespace shimmer_teqp
{
  namespace utilities
  {
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
    template<typename mol_fracs_type,
             typename c_vec_type,
             typename value_type>
    auto pseudo_critical_point(const mol_fracs_type& mol_fracs,
                               const c_vec_type& Tc_vec,
                               const c_vec_type& vc_vec,
                               const value_type& tolerance)
    {
      using size_type = typename c_vec_type::size_type;

      assert(static_cast<size_type>(mol_fracs.size()) == Tc_vec.size() &&
             vc_vec.size() == Tc_vec.size());

      gerg_data::Pseudo_critical_point<value_type> pseudo_critical_point;

      pseudo_critical_point.Tcx = value_type();
      pseudo_critical_point.Vcx = value_type();
      pseudo_critical_point.Dcx = value_type();

      for (size_type i = 0; i < Tc_vec.size(); ++i)
      {
        pseudo_critical_point.Tcx += mol_fracs[i] *
                                     Tc_vec[i];
        pseudo_critical_point.Vcx += mol_fracs[i] *
                                     vc_vec[i];
      }

      if (pseudo_critical_point.Vcx > tolerance)
        pseudo_critical_point.Dcx = 1.0 / pseudo_critical_point.Vcx;

      return pseudo_critical_point;
    }
    // *********************************************************
  }
}

#endif
