#pragma once

#include <iostream>
#include <vector>

namespace shimmer {

struct setting_remi_wo {
    int     u_snum;     // User number of the station
    int     i_snum;     // Internal number of the station
    double  Lmin;       // Minimum mass flow rate
    double  Lmax;       // Maximum mass flow rate
    double  Pmin;       // Minumum pressure
    double  Pmax;       // Maximum pressure

    std::vector<sample> Pprofile;

    bool operator<(const setting_remi_wo& other) {
        return i_snum < other.i_snum;
    }
};

inline std::ostream&
operator<<(std::ostream& os, const setting_remi_wo& s) {
    os << "ReMi w/o backflow - unum: " << s.u_snum << ", inum: " << s.i_snum << " - ";
    os << "Lmin = " << s.Lmin << ", Lmax = " << s.Lmax << ", Pmin = " << s.Pmin;
    os << ", Pmax = " << s.Pmax << ", profile samples: " << s.Pprofile.size();
    return os;
}

} // namespace shimmer