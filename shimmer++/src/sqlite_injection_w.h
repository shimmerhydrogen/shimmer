#pragma once

#include <iostream>
#include <vector>
#include "sqlite.hpp"

namespace shimmer {

struct setting_injection_w {
    int     u_snum;     // User number of the station
    int     i_snum;     // Internal number of the station
    double  Lmin;       // Minimum mass flow rate
    double  Lmax;       // Maximum mass flow rate
    double  Pmin;       // Minumum pressure
    double  Pmax;       // Maximum pressure
    double  f;          // scaling parameter in [0,1]

    std::vector<sample> Pprofile;

    bool operator<(const setting_injection_w& other) {
        return i_snum < other.i_snum;
    }
};

inline std::ostream&
operator<<(std::ostream& os, const setting_injection_w& s) {
    os << "ReMi w/o backflow - unum: " << s.u_snum << ", inum: " << s.i_snum << " - ";
    os << "Lmin = " << s.Lmin << ", Lmax = " << s.Lmax << ", Pmin = " << s.Pmin;
    os << ", Pmax = " << s.Pmax << ", scale = " << s.f << ", profile samples: ";
    os << s.Pprofile.size();
    return os;
}

} // namespace shimmer