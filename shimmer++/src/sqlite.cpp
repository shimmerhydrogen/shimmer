#include <iostream>
#include <cstdio>

#include "sqlite.hpp"

namespace shimmer {

network_db::network_db()
    : db_(nullptr)
{}

network_db::network_db(const std::string& filename)
{
    open(filename);
}

int
network_db::open(const std::string& filename)
{
    int rc;
    rc = sqlite3_open(filename.c_str(), &db_);
    if(rc) {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_));
        sqlite3_close(db_);
        return 1;
    }

    return 0;
}

network_db::~network_db()
{
    sqlite3_close(db_);
}

int
network_db::read_stations()
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM stations";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        //int num_cols = sqlite3_column_count(stmt);
        vertex_properties vp;
        vp.name = (char *) sqlite3_column_text(stmt, 0);
        vp.node_num = sqlite3_column_int(stmt, 1);
        vp.height = sqlite3_column_double(stmt, 2);
        vmap_[vp.node_num] = std::move(vp);
    }

    std::cout << "Read " << vmap_.size() << " stations" << std::endl;

    sqlite3_finalize(stmt);

    return 0;
}

int
network_db::read_station_fd_parameters()
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM station_fd_parameters";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        //int num_cols = sqlite3_column_count(stmt);
        int node_num = sqlite3_column_int(stmt, 0);
        vmap_[node_num].pressure = sqlite3_column_double(stmt, 1);
        vmap_[node_num].mass_flow = sqlite3_column_double(stmt, 2);
    }

    std::cout << "Read " << vmap_.size() << " stations" << std::endl;

    sqlite3_finalize(stmt);

    for (auto& [n, v] : vmap_)
        std::cout << v << std::endl;

    return 0;
}

int
network_db::read_pipelines()
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM pipelines";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    size_t i = 0;
    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int from = sqlite3_column_int(stmt, 1);
        int to = sqlite3_column_int(stmt, 2);
        int type = sqlite3_column_int(stmt, 3);
        i++;
    }

    std::cout << "Read " << i << " pipelines" << std::endl;

    sqlite3_finalize(stmt);

    return 0;
}

int
network_db::populate(infrastructure_graph& g)
{
    if (!db_)
        return 1;

    // check errors!!!!!
    read_stations();
    read_station_fd_parameters();
    read_pipelines();

    return 0;
}

} // namespace shimmer