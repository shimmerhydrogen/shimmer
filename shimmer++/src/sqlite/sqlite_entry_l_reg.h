#pragma once

#include <iostream>
#include <vector>

namespace shimmer {

struct setting_entry_l_reg {
    int     u_snum;     // User number of the station
    int     i_snum;     // Internal number of the station
    double  Lmin;       // Minimum mass flow rate
    double  Lmax;       // Maximum mass flow rate
    double  Pmin;       // Minumum pressure
    double  Pmax;       // Maximum pressure
    double  f;          // scaling parameter in [0,1]

    std::vector<sample> Pprofile;
    std::vector<sample> Lprofile;

    bool operator<(const setting_entry_l_reg& other) {
        return i_snum < other.i_snum;
    }
};

inline std::ostream&
operator<<(std::ostream& os, const setting_entry_l_reg& s) {
    os << "Injection w/ press. control - unum: " << s.u_snum << ", inum: " << s.i_snum << " - ";
    os << "Lmin = " << s.Lmin << ", Lmax = " << s.Lmax << ", Pmin = " << s.Pmin;
    os << ", Pmax = " << s.Pmax << ", scale = " << s.f << ", profile samples: ";
    os << s.Pprofile.size();
    return os;
}

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_entry_l_reg>& settings);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_entry_l_reg>& settings);

} //namespace database

} // namespace shimmer