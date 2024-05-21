#pragma once

#include <sqlite3.h>
#include "infrastructure_graph.h"

namespace shimmer {

int db_read_stations(sqlite3 *, infrastructure_graph&);

class network_db {

    sqlite3 *db_;
    bool verbose_;
    std::map<int, vertex_properties> vmap_;

    int read_stations();
    int read_station_fd_parameters();

public:
    network_db();
    network_db(const std::string&);
    ~network_db();

    int open(const std::string&);
    int populate(infrastructure_graph&);
};

} //namespace shimmer