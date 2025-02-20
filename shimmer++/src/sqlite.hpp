#pragma once

#include <optional>
#include <sqlite3.h>
#include "infrastructure_graph.h"

namespace shimmer {
struct sample {
    double  time;
    double  value;

    bool operator<(const sample& other) {
        return time < other.time;
    }
};

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
#include "sqlite_conspoint_wo.h"

namespace shimmer {

class network_database {

    sqlite3 *db_;
    bool verbose_;

    /* BEGIN I have the impression that this stuff does not belong here */
    std::vector<std::optional<int>>                 s_u2i;
    std::vector<int>                                s_i2u;
    std::vector<std::optional<vertex_descriptor>>   s_u2vd;
    std::vector<vertex_descriptor>                  s_i2vd;

    std::vector<setting_outlet>         settings_outlet;
    std::vector<setting_entry_p_reg>        settings_entry_p_reg;
    std::vector<setting_entry_l_reg>    settings_entry_l_reg;
    std::vector<setting_conspoint_wo>   settings_conspoint_wo;

    /* END I have the impression that this stuff does not belong here */

    int import_stations(infrastructure_graph& g);
    int import_station_parameters(infrastructure_graph& g);
    int import_pipelines(infrastructure_graph& g);
    int renumber_stations();


    int import_outlet(std::vector<setting_outlet>&);
    int import_entry_p_reg(std::vector<setting_entry_p_reg>&);
    int import_entry_l_reg(std::vector<setting_entry_l_reg>&);
    int import_conspoint_wo(std::vector<setting_conspoint_wo>&);

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
