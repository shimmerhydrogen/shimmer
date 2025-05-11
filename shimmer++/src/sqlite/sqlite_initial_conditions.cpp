#include <optional>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

namespace shimmer {

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<station_initial_condition>& sics)
{
    int rc;
    std::string q = "SELECT * FROM station_initial_conditions";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int u_snum = sqlite3_column_int(stmt, 0);
        int i_snum = convert_u2i(s_u2i, u_snum).value();
        station_initial_condition sic;
        sic.i_snum = i_snum;
        sic.init_P = sqlite3_column_double(stmt, 1);
        sic.init_L = sqlite3_column_double(stmt, 2);
        sics.push_back(sic);
    }

    std::sort(sics.begin(), sics.end());

    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;
}

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<station_initial_condition>& sics)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO station_initial_conditions "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (int i = 0; i < sics.size(); i++) {
        rc = sqlite3_bind_int(stmt, 1, convert_i2u(s_i2u, sics[i].i_snum));
        rc = sqlite3_bind_double(stmt, 2, sics[i].init_P);
        rc = sqlite3_bind_double(stmt, 3, sics[i].init_L);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<pipe_initial_condition>& pics)
{
    int rc;
    std::string q = "SELECT * FROM pipe_initial_conditions";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int u_sfrom = sqlite3_column_int(stmt, 1);
        int i_sfrom = convert_u2i(s_u2i, u_sfrom).value();
        int u_sto = sqlite3_column_int(stmt, 2);
        int i_sto = convert_u2i(s_u2i, u_sto).value();

        pipe_initial_condition pic;
        pic.i_sfrom = i_sfrom;
        pic.i_sto = i_sto;
        pic.init_G = sqlite3_column_double(stmt, 3);
        pics.push_back(pic);
    }

    std::sort(pics.begin(), pics.end());

    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;
}

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<pipe_initial_condition>& pics)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO pipe_initial_conditions "
        "VALUES (?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (int i = 0; i < pics.size(); i++) {
        auto from = convert_i2u(s_i2u, pics[i].i_sfrom);
        auto to = convert_i2u(s_i2u, pics[i].i_sto);
        auto init_G = pics[i].init_G;
        std::string pipe_name = "pipe-" + std::to_string(from) + "-" + std::to_string(to);

        rc = sqlite3_bind_text(stmt, 1, pipe_name.c_str(), pipe_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, from);
        rc = sqlite3_bind_int(stmt, 3, to);
        rc = sqlite3_bind_int(stmt, 4, init_G);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    return 0;
}

} // namespace database

} //namespace shimmer