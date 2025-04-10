#pragma once

#include <iostream>
#include <vector>
#include "sqlite.hpp"

namespace shimmer {


struct setting_pipe {
    int         i_sfrom;
    int         i_sto;

    double      diameter;
    double      length;
    double      roughness;

    bool operator<(const setting_pipe& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

} // namespace shimmer