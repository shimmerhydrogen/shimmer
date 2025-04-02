/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 *
 * The authors (C) 2023, 2024, 2025
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
    auto ppitor = std::lower_bound(ppsbegin, ppsend, {i_sfrom, i_sto}, comp);
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

} // namespace shimmer

#include "sqlite_outlet.h"
#include "sqlite_entry_p_reg.h"
#include "sqlite_entry_l_reg.h"
#include "sqlite_exit_l_reg.h"
#include "sqlite_initial_conditions.h"
#include "sqlite_compressor.h"

namespace shimmer {

template<typename T>
using optvector = std::vector<std::optional<T>>;

class network_database {

    sqlite3 *db_;
    bool verbose_;

    /* BEGIN I have the impression that this stuff does not belong here */
    optvector<int>                      s_u2i;
    std::vector<int>                    s_i2u;
    optvector<vertex_descriptor>        s_u2vd;
    std::vector<vertex_descriptor>      s_i2vd;

    std::vector<setting_outlet>         settings_outlet;
    std::vector<setting_entry_p_reg>    settings_entry_p_reg;
    std::vector<setting_entry_l_reg>    settings_entry_l_reg;
    std::vector<setting_exit_l_reg>     settings_exit_l_reg;
    std::vector<setting_compr_stat>     settings_compr_stat;

    std::vector<station_initial_condition>  sics;
    std::vector<pipe_initial_condition>     pics;
    /* END I have the impression that this stuff does not belong here */

    int import_stations(infrastructure_graph&);
    int import_station_parameters(infrastructure_graph&);
    int import_pipelines(infrastructure_graph&);
    int renumber_stations();

    int import_outlet(std::vector<setting_outlet>&);
    int import_entry_p_reg(std::vector<setting_entry_p_reg>&);
    int import_entry_l_reg(std::vector<setting_entry_l_reg>&);
    int import_exit_l_reg(std::vector<setting_exit_l_reg>&);

    int import_compr_stat(std::vector<setting_compr_stat>&);

    int import_station_initial_conditions();
    int import_pipe_initial_conditions();

    int populate_type_dependent_station_data(vertex_properties&);

    std::optional<table_name_pair_t> limits_and_profile_table_names(station_type);

public:
    network_database();
    network_database(const std::string&);
    ~network_database();

    int open(const std::string&);
    int populate_graph(infrastructure_graph&);
};

} //namespace shimmer
