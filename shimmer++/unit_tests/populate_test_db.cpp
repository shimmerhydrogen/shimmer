#include <iostream>
#include <sstream>
#include <vector>
#include <tuple>
#include <optional>

#include <sqlite3.h>

#include "sqlite/sqlite.hpp"
#include "network_elements.h"

namespace shimmer::testing
{

struct state
{
    sqlite3 *db;
};

int populate_stations(const state& st, const std::vector<station_type>& stations)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO stations(s_number, s_name, t_type) "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < stations.size(); i++) {
        std::string station_name = "station" + std::to_string(i);

        rc = sqlite3_bind_int(stmt, 1, i);
        rc = sqlite3_bind_text(stmt, 2, station_name.c_str(), station_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 3, +stations[i]);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);
    
    return 0;
}

using edge_type = std::tuple<int, int, shimmer::pipe_type>;
int populate_pipes(const state& st, const std::vector<edge_type>& pipes)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO pipelines(p_name, s_from, s_to, p_type) "
        "VALUES (?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < pipes.size(); i++) {
        auto from = std::get<0>(pipes[i]);
        auto to = std::get<1>(pipes[i]);
        auto type = std::get<2>(pipes[i]);
        std::string pipe_name = "pipe-" + std::to_string(from) + "-" + std::to_string(to);

        rc = sqlite3_bind_text(stmt, 1, pipe_name.c_str(), pipe_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, from);
        rc = sqlite3_bind_int(stmt, 3, to);
        rc = sqlite3_bind_int(stmt, 4, +type);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

struct plain_pipe_params {
    int     from;
    int     to;
    double  length;
    double  diameter;
    double  friction_factor;
};

int populate_plain_pipe_params(const state& st, const std::vector<plain_pipe_params>& pp)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO pipe_parameters "
        "VALUES (?, ?, ?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < pp.size(); i++) {
        auto from = pp[i].from;
        auto to = pp[i].to;
        std::string pipe_name = "pipe-" + std::to_string(from) + "-" + std::to_string(to);

        rc = sqlite3_bind_text(stmt, 1, pipe_name.c_str(), pipe_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, from);
        rc = sqlite3_bind_int(stmt, 3, to);
        rc = sqlite3_bind_double(stmt, 4, pp[i].diameter);
        rc = sqlite3_bind_double(stmt, 5, pp[i].length);
        rc = sqlite3_bind_double(stmt, 6, pp[i].friction_factor);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

using stat_ic_type = std::tuple<int, double, double>;
int populate_station_ics(const state& st, const std::vector<stat_ic_type>& sics)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO station_initial_conditions "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < sics.size(); i++) {
        auto s_number = std::get<0>(sics[i]);
        auto init_P = std::get<1>(sics[i]);
        auto init_L = std::get<2>(sics[i]);

        rc = sqlite3_bind_int(stmt, 1, s_number);
        rc = sqlite3_bind_double(stmt, 2, init_P);
        rc = sqlite3_bind_double(stmt, 3, init_L);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

using pipe_ic_type = std::tuple<int, int, double>;
int populate_pipe_ics(const state& st, const std::vector<pipe_ic_type>& pipe_ics)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO pipe_initial_conditions "
        "VALUES (?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < pipe_ics.size(); i++) {
        auto from = std::get<0>(pipe_ics[i]);
        auto to = std::get<1>(pipe_ics[i]);
        auto init_G = std::get<2>(pipe_ics[i]);
        std::string pipe_name = "pipe-" + std::to_string(from) + "-" + std::to_string(to);

        rc = sqlite3_bind_text(stmt, 1, pipe_name.c_str(), pipe_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, from);
        rc = sqlite3_bind_int(stmt, 3, to);
        rc = sqlite3_bind_int(stmt, 4, init_G);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int populate_limits(const state& st, int station, station_type type)
{
    auto names_opt = limits_and_profile_table_names(st.db, type);
 
    if ( not names_opt ) {
        std::cerr << "    Not storing limits data for station type " << +type << std::endl;
        return -1;
    }

    auto [limits_tab, profiles_tab] = names_opt.value();

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + limits_tab + " (s_number, lim_Lmin, lim_Lmax, lim_Pmin, lim_Pmax) "
        "VALUES (?, ?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    double lim_Pmin = 60e5;
    double lim_Pmax = 80e5;
    double lim_Lmin = -300;
    double lim_Lmax = -10;

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

        rc = sqlite3_bind_int(stmt, 1, station);
        rc = sqlite3_bind_double(stmt, 2, lim_Lmin);
        rc = sqlite3_bind_double(stmt, 3, lim_Lmax);
        rc = sqlite3_bind_double(stmt, 4, lim_Pmin);
        rc = sqlite3_bind_double(stmt, 5, lim_Pmax);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
     
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int populate_entry_l_reg_profile(const state& st, int station,
    const std::vector<double>& P, const std::vector<double>& L,
    double delta_t)
{
    auto type = station_type::ENTRY_L_REG;
    std::cout << "  L profile for type " << +type << " station " << station << std::endl;
    auto names_opt = limits_and_profile_table_names(st.db, type);
 
    if ( not names_opt ) {
        std::cerr << "    Not storing profile data for station type " << +type << std::endl;
        return -1;
    }

    auto [limits_tab, profiles_tab] = names_opt.value();

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + profiles_tab + " (s_number, prf_time, prf_Pset, prf_Lset) "
        "VALUES (?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < L.size(); i++) {
        rc = sqlite3_bind_int(stmt, 1, station);
        rc = sqlite3_bind_int(stmt, 2, i*delta_t);
        rc = sqlite3_bind_int(stmt, 3, P[i]);
        rc = sqlite3_bind_int(stmt, 4, L[i]);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int populate_P_profile(const state& st, int station, station_type type,
    const std::vector<double>& P, double delta_t)
{
    std::cout << "  L profile for type " << +type << " station " << station << std::endl;
    auto names_opt = limits_and_profile_table_names(st.db, type);
 
    if ( not names_opt ) {
        std::cerr << "    Not storing profile data for station type " << +type << std::endl;
        return -1;
    }

    auto [limits_tab, profiles_tab] = names_opt.value();

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + profiles_tab + " (s_number, prf_time, prf_Pset) "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < P.size(); i++) {
        rc = sqlite3_bind_int(stmt, 1, station);
        rc = sqlite3_bind_int(stmt, 2, i*delta_t);
        rc = sqlite3_bind_int(stmt, 3, P[i]);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int populate_L_profile(const state& st, int station, station_type type,
    const std::vector<double>& L, double delta_t)
{
    std::cout << "  L profile for type " << +type << " station " << station << std::endl;
    auto names_opt = limits_and_profile_table_names(st.db, type);
 
    if ( not names_opt ) {
        std::cerr << "    Not storing profile data for station type " << +type << std::endl;
        return -1;
    }

    auto [limits_tab, profiles_tab] = names_opt.value();

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + profiles_tab + " (s_number, prf_time, prf_Lset) "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(st.db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(st.db, "BEGIN TRANSACTION", 0, 0, 0);

    for (int i = 0; i < L.size(); i++) {
        rc = sqlite3_bind_int(stmt, 1, station);
        rc = sqlite3_bind_int(stmt, 2, i*delta_t);
        rc = sqlite3_bind_int(stmt, 3, L[i]);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(st.db, "COMMIT", 0, 0, 0);

    rc = sqlite3_finalize(stmt);

    return 0;
}

int initialize_test_db(void)
{
    state st;
    int rc = sqlite3_open("test_network_1.db", &st.db);
    if (rc) {
        std::cerr << "Can't open database: ";
        std::cerr << sqlite3_errmsg(st.db) << std::endl;
        return -1;
    }

    /* Stations */
    std::cout << "Stations" << std::endl;
    std::vector<station_type> station_types {
        station_type::ENTRY_P_REG,
        station_type::PRIVATE_OUTLET,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::EXIT_L_REG,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::EXIT_L_REG,
        station_type::JUNCTION,
        station_type::EXIT_L_REG,
        station_type::ENTRY_L_REG,
        station_type::PRIVATE_OUTLET
    };
    populate_stations(st, station_types);
    
    /* Pipelines */
    std::cout << "Pipelines" << std::endl;
    std::vector<edge_type> pipes = {
        {  0,  3, pipe_type::PIPE },
        {  1,  2, pipe_type::PIPE },
        {  2,  3, pipe_type::PIPE },
        {  2,  4, pipe_type::PIPE },
        {  3,  4, pipe_type::PIPE },
        {  4,  5, pipe_type::PIPE },
        {  4,  7, pipe_type::PIPE },
        {  6,  4, pipe_type::PIPE },
        {  7,  6, pipe_type::PIPE },
        { 11,  6, pipe_type::PIPE },
        { 12,  7, pipe_type::PIPE },
        {  8,  7, pipe_type::PIPE },
        {  7,  9, pipe_type::PIPE },
        {  9, 10, pipe_type::PIPE },
        {  3,  9, pipe_type::PIPE }
    };
    populate_pipes(st, pipes);

    std::vector<plain_pipe_params> pp = {
        {  0,  3,   80000,	1.2,	1.20E-05 },
        {  1,  2,   16000,	0.6,	1.20E-05 },
        {  2,  3,   40000,	0.8,	1.20E-05 },
        {  2,  4,  160000,	0.7,	1.20E-05 },
        {  3,  4,  200000,	0.8,	1.20E-05 },
        {  4,  5,   24000,	0.6,	1.20E-05 },
        {  4,  7,  120000,	0.2,	1.20E-05 },
        {  6,  4,   80000,	0.9,	1.20E-05 },
        {  7,  6,   64000,	0.7,	1.20E-05 },
        { 11,  6,  240000,	0.6,	1.20E-05 },
        { 12,  7,   28000,	0.2,	1.20E-05 },
        {  8,  7,   80000,	0.9,	1.20E-05 },
        {  7,  9,  160000,	0.7,	1.20E-05 },
        {  9, 10,   40000,	0.3,	1.20E-05 },
        {  3,  9,  320000,	0.9,	1.20E-05 }
    };
    populate_plain_pipe_params(st, pp);


    /* Initial conditions for stations */
    std::cout << "Initial conditions for stations" << std::endl;
    std::vector<stat_ic_type> stat_ics = {
        {  0, 70.000000000000000*1e5, -2.300000000000000e+02 },
        {  1, 70.000000000000000*1e5, -2.280833333333333e+02 },
        {  2, 69.299999999999997*1e5, -2.261666666666667e+02 },
        {  3, 69.299999999999997*1e5, -2.242500000000000e+02 },
        {  4, 68.606999999999999*1e5, -2.223333333333334e+02 },
        {  5, 67.920929999999998*1e5, -2.204166666666667e+02 },
        {  6, 67.241720700000002*1e5, -2.185000000000000e+02 },
        {  7, 67.920929999999998*1e5, -2.204166666666667e+02 },
        {  8, 67.241720700000002*1e5, -2.223333333333333e+02 },
        {  9, 67.241720700000002*1e5, -2.242500000000000e+02 },
        { 10, 66.569303493000007*1e5, -2.261666666666666e+02 },
        { 11, 70.000000000000000*1e5, -2.280833333333333e+02 },
        { 12, 67.241720700000002*1e5, -2.300000000000000e+02 }
    };
    populate_station_ics(st, stat_ics);

    /* Initial conditions for pipes */
    std::cout << "Initial conditions for pipes" << std::endl;
    std::vector<pipe_ic_type> pipe_ics;
    pipe_ics.reserve(pipes.size());
    for (int i = 0; i < pipes.size(); i++)
        pipe_ics.push_back( {std::get<0>(pipes[i]), std::get<1>(pipes[i]), 50} );
    populate_pipe_ics(st, pipe_ics);

    /* Limits */
    std::cout << "Limits" << std::endl;
    for (int i = 0; i < station_types.size(); i++)
        populate_limits(st, i, station_types[i]);

    /* Profiles */
    std::cout << "Profiles" << std::endl;
    std::vector<double> node0_P = {
        -75.000000000000000, -74.375000000000000, -73.750000000000000,
        -73.125000000000000, -72.500000000000014, -71.875000000000014,
        -71.250000000000000, -71.875000000000000, -72.499999999999986,
        -73.124999999999986, -73.749999999999986, -74.374999999999986,
        -75.000000000000000, -77.499999999999986, -79.999999999999986,
        -82.499999999999972, -84.999999999999972, -87.499999999999957,
        -90.000000000000000, -87.499999999999986, -85.000000000000000,
        -82.500000000000000, -80.000000000000014, -77.500000000000028,
        -7.5000000e+01 };
    populate_P_profile(st, 0, station_types[0], node0_P, 3600);

    std::vector<double> node1_G = {
        20.000000000000000,  19.166666666666668,  18.333333333333336,
        17.500000000000004,  16.666666666666671,  15.833333333333337,
        15.000000000000000,  16.666666666666668,  18.333333333333336,
        20.000000000000000,  21.666666666666664,  23.333333333333329,
        25.000000000000000,  23.333333333333336,  21.666666666666671,
        20.000000000000004,  18.333333333333336,  16.666666666666671,
        15.000000000000000,  16.666666666666668,  18.333333333333336,
        20.000000000000000,  21.666666666666664,  23.333333333333329,
        2.5000000e+01 };
    populate_L_profile(st, 1, station_types[1], node1_G, 3600);

    std::vector<double> node5_G = {
        20.000000000000000,  19.333333333333332,  18.666666666666668,
        18.000000000000000,  17.333333333333336,  16.666666666666668,
        16.000000000000000,  16.666666666666668,  17.333333333333336,
        18.000000000000000,  18.666666666666668,  19.333333333333332,
        20.000000000000000,  19.333333333333332,  18.666666666666668,
        18.000000000000000,  17.333333333333336,  16.666666666666668,
        16.000000000000000,  16.666666666666668,  17.333333333333336,
        18.000000000000000,  18.666666666666668,  19.333333333333332,
        2.0000000e+01 };
    populate_L_profile(st, 5, station_types[5], node5_G, 3600);

    std::vector<double> node8_G = {
        50.000000000000000,  52.083333333333336,  54.166666666666671,
        56.250000000000014,  58.333333333333350,  60.416666666666686,
        62.500000000000000,  58.750000000000000,  55.000000000000007,
        51.250000000000007,  47.500000000000007,  43.750000000000014,
        40.000000000000000,  40.000000000000000,  40.000000000000000,
        40.000000000000000,  40.000000000000000,  40.000000000000000,
        40.000000000000000,  41.666666666666671,  43.333333333333336,
        45.000000000000000,  46.666666666666664,  48.333333333333336,
        5.0000000e+01 };
    populate_L_profile(st, 8, station_types[8], node8_G, 3600);

    std::vector<double> node10_G = {
        15.000000000000000,  14.750000000000000,  14.500000000000002,
        14.250000000000004,  14.000000000000004,  13.750000000000004,
        13.500000000000000,  14.250000000000000,  15.000000000000000,
        15.750000000000000,  16.500000000000000,  17.250000000000004,
        18.000000000000000,  17.499999999999996,  17.000000000000000,
        16.500000000000000,  16.000000000000004,  15.500000000000005,
        15.000000000000000,  15.250000000000000,  15.499999999999998,
        15.749999999999996,  15.999999999999996,  16.249999999999996,
        1.6500000e+01 } ;
    populate_L_profile(st, 10, station_types[10], node10_G, 3600);

    std::vector<double> node11_P(25, 7500000.0);
    std::vector<double> node11_L = {
        -40.000000000000000, -38.666666666666664, -37.333333333333336,
        -36.000000000000000, -34.666666666666671, -33.333333333333336,
        -32.000000000000000, -33.333333333333336, -34.666666666666671,
        -36.000000000000000, -37.333333333333336, -38.666666666666664,
        -40.000000000000000, -41.333333333333329, -42.666666666666657,
        -43.999999999999986, -45.333333333333314, -46.666666666666643,
        -48.000000000000000, -48.333333333333329, -48.666666666666664,
        -48.999999999999993, -49.333333333333329, -49.666666666666657,
        -5.0000000e+01 };
    populate_entry_l_reg_profile(st, 11, node11_P, node11_L, 3600);

    std::vector<double> node12_G = {
        10.000000000000000,   9.666666666666666,   9.333333333333334,
         9.000000000000000,   8.666666666666668,   8.333333333333334,
         8.000000000000000,   8.333333333333334,   8.666666666666668,
         9.000000000000000,   9.333333333333334,   9.666666666666666,
        10.000000000000000,   9.833333333333334,   9.666666666666668,
         9.500000000000002,   9.333333333333336,   9.166666666666670,
         9.000000000000000,   9.583333333333334,  10.166666666666666,
        10.750000000000000,  11.333333333333332,  11.916666666666666,
         1.2500000e+01 };
    populate_L_profile(st, 12, station_types[12], node12_G, 3600);

    std::cout << "Closing DB" << std::endl;
    sqlite3_close(st.db);

    return 0;
}

} //namespace shimmer::testing



    

int main(void)
{
    shimmer::testing::initialize_test_db();
    return 0;
}