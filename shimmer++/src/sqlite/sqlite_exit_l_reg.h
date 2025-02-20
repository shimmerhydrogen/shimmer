#pragma once

#include <iostream>
#include <vector>

namespace shimmer {

struct setting_exit_l_reg {
    int     u_snum;     // User number of the station
    int     i_snum;     // Internal number of the station
    double  Lmin;       // Minimum mass flow rate
    double  Lmax;       // Maximum mass flow rate
    double  Pmin;       // Minumum pressure
    double  Pmax;       // Maximum pressure

    std::vector<sample> Lprofile;

    bool operator<(const setting_exit_l_reg& other) {
        return i_snum < other.i_snum;
    }
};

inline std::ostream&
operator<<(std::ostream& os, const setting_exit_l_reg& s) {
    os << "Consumption point w/o P control - unum: " << s.u_snum << ", inum: " << s.i_snum << " - ";
    os << "Lmin = " << s.Lmin << ", Lmax = " << s.Lmax << ", Pmin = " << s.Pmin;
    os << ", Pmax = " << s.Pmax << ", profile samples: " << s.Lprofile.size();
    return os;
}

} // namespace shimmer