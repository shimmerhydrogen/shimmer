/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <optional>
#include <sqlite3.h>
#include "infrastructure_graph.h"

namespace shimmer {

template<typename T>
concept station_props = requires {
    T::i_snum;
};

template<station_props S>
auto
lookup(const std::vector<S>& sps, int i_snum)
{
    auto comp = [](const S& sp, int i) {
        return sp.i_snum < i;
    };

    auto spsbegin = sps.begin();
    auto spsend = sps.end();
    auto spitor = std::lower_bound(spsbegin, spsend, i_snum, comp);
    if ( (spitor == spsend) or ((*spitor).i_snum != i_snum) )
        return spsend;

    return spitor;
}

template<typename T>
concept pipe_props = requires {
    T::i_sfrom;
    T::i_sto;
};

template<pipe_props P>
auto
lookup(const std::vector<P>& pps, int i_sfrom, int i_sto)
{
    auto comp = [](const P& pp, const std::pair<int, int>& fromto) {
        return std::pair{pp.i_sfrom, pp.i_sto} < fromto;
    };

    auto ppsbegin = pps.begin();
    auto ppsend = pps.end();
    auto ppitor = std::lower_bound(ppsbegin, ppsend, std::pair{i_sfrom, i_sto}, comp);
    if ( (ppitor == ppsend) or ((*ppitor).i_sfrom != i_sfrom)
        or ((*ppitor).i_sto != i_sto) )
        return ppsend;

    return ppitor;
}
struct sample {
    double  time;
    double  value;

    bool operator<(const sample& other) {
        return time < other.time;
    }
};

inline double
interp(const std::vector<sample>& samples, double time) {
    if (samples.size() == 0)
        return 0.0;

    const auto& first = samples.front();
    if (time < first.time)
        return first.value;

    for (size_t i = 1; i < samples.size(); i++) {
        const auto& t0 = samples[i-1].time;
        const auto& v0 = samples[i-1].value;
        const auto& t1 = samples[i].time;
        const auto& v1 = samples[i].value;

        if ( (time >= t0) and (time < t1) ) {
            return v0 + (time - t0)*(v1-v0)/(t1 - t0);
        }
    }
    
    const auto& last = samples.back();
    return last.value;
}

inline std::ostream&
operator<<(std::ostream& os, const sample& s) {
    os << "(" << s.time << ", " << s.value << ")";
    return os;
}

using table_name_pair_t = std::pair<std::string, std::string>;
std::optional<table_name_pair_t>
limits_and_profile_table_names(sqlite3 *, station_type);

enum class setting_table {
    limits,
    profiles
};
std::optional<std::string>
table_name(sqlite3 *db, setting_table st, station_type stat_type);

template<typename T>
using optvector = std::vector<std::optional<T>>;

inline int
convert_i2u(const std::vector<int>& s_i2u, int idx)
{
    if (s_i2u.size() > 0)
        return s_i2u.at(idx);

    return idx;
}

inline std::optional<int>
convert_u2i(const optvector<int>& s_u2i, int idx)
{
    if (s_u2i.size() > 0 and idx < s_u2i.size())
        return s_u2i[idx];

    return idx;
}

} // namespace shimmer

#include "sqlite_outlet.h"
#include "sqlite_entry_p_reg.h"
#include "sqlite_entry_l_reg.h"
#include "sqlite_exit_l_reg.h"
#include "sqlite_initial_conditions.h"
#include "sqlite_pipe.h"
#include "sqlite_compressor.h"
#include "sqlite_reduction.h"
#include "sqlite_gases.h"
#include "sqlite_infra.h"


namespace shimmer {

template<typename T>
concept limits = requires {
    T::Pmin;
    T::Pmax;
    T::Lmin;
    T::Lmax;
};

/* Temporary helper */
auto
convert_limits(const limits auto& l)
{
    std::vector<pair_input_t> user_constraints = {
        std::make_pair(P_GREATER_EQUAL, l.Pmin),
        std::make_pair(P_LOWER_EQUAL,   l.Pmax),
        std::make_pair(L_GREATER_EQUAL, l.Lmin),
        std::make_pair(L_LOWER_EQUAL,   l.Lmax)
    };

    return user_constraints;
}

template<typename T>
concept profile =
    requires { T::Lprofile; } ||
    requires { T::Pprofile; };

/* Temporary helper */
auto
convert_Lprof(const profile auto& p)
{
    Eigen::VectorXd ret( p.Lprofile.size() );
    for (size_t i = 0; i < p.Lprofile.size(); i++)
        ret[i] = p.Lprofile[i].value;
    return ret;
}

/* Temporary helper */
auto
convert_Pprof(const profile auto& p)
{
    Eigen::VectorXd ret( p.Pprofile.size() );
    for (size_t i = 0; i < p.Pprofile.size(); i++)
        ret[i] = p.Pprofile[i].value;
    return ret;
}

} //namespace shimmer
