#pragma once

#include <iostream>
#include <vector>

namespace shimmer {

struct setting_outlet {
    int     i_snum;     // Internal number of the station

    std::vector<sample> Lprofile;

    bool operator<(const setting_outlet& other) {
        return i_snum < other.i_snum;
    }
};

inline std::ostream&
operator<<(std::ostream& os, const setting_outlet& s) {
    os << "Outlet - inum: " << s.i_snum << " - ";
    os << "profile samples: ";
    os << s.Lprofile.size();
    return os;
}

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_outlet>& settings);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_outlet>& settings);

} //namespace database

} // namespace shimmer