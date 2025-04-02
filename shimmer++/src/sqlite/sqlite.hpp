#pragma once

#include <optional>
#include <sqlite3.h>
#include "infrastructure_graph.h"

namespace shimmer {

struct station_initial_condition {
    int     s_number;
    double  init_P;
    double  init_L;

    bool operator<(const station_initial_condition& other) const {
        return s_number < other.s_number;
    }
};

struct pipe_initial_condition {
    int     s_from;
    int     s_to;
    double  init_G;

    bool operator<(const pipe_initial_condition& other) const {
        return (s_from < other.s_from) or
            (s_from == other.s_from and s_to < other.s_to);
    }
};

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

    const auto& last = samples.back();
    if (time >= last.time)
        return last.value;

    for (size_t i = 1; i < samples.size(); i++) {
        const auto& t0 = samples[i-1].time;
        const auto& v0 = samples[i-1].value;
        const auto& t1 = samples[i].time;
        const auto& v1 = samples[i].value;

        if ( (time >= t0) and (time < t1) ) {
            return v0 + (time - t0)*(v1-v0)/(t1 - t0);
        }
    }
    
    __builtin_unreachable();
}

inline std::ostream&
operator<<(std::ostream& os, const sample& s) {
    os << "(" << s.time << ", " << s.value << ")";
    return os;
}

using table_name_pair_t = std::pair<std::string, std::string>;

} // namespace shimmer

template<typename T>
    requires std::is_enum_v<T>
constexpr auto operator+(T e) {
    return std::underlying_type_t<T>(e);
}

#include "sqlite_outlet.h"
#include "sqlite_entry_p_reg.h"
#include "sqlite_entry_l_reg.h"
#include "sqlite_exit_l_reg.h"

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

    int import_station_initial_conditions();
    int import_pipe_initial_conditions();

    int populate_type_dependent_station_data(vertex_properties&);

    std::optional<table_name_pair_t> limits_and_profile_table_names(int stat_type);

public:
    network_database();
    network_database(const std::string&);
    ~network_database();

    int open(const std::string&);
    int populate_graph(infrastructure_graph&);
};

} //namespace shimmer
