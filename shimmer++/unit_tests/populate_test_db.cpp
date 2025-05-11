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

    std::vector<setting_pipe> settings_pipe = {
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
    database::store(st.db, {}, settings_pipe);



    /* Initial conditions for stations */
    std::cout << "Initial conditions for stations" << std::endl;
    std::vector<station_initial_condition> stat_ics = {
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
    database::store(st.db, {}, stat_ics);

    /* Initial conditions for pipes */
    std::cout << "Initial conditions for pipes" << std::endl;
    std::vector<pipe_initial_condition> pipe_ics;
    pipe_ics.reserve(pipes.size());
    for (int i = 0; i < pipes.size(); i++)
        pipe_ics.push_back( {std::get<0>(pipes[i]), std::get<1>(pipes[i]), 50} );
    database::store(st.db, {}, pipe_ics);

    /* Profiles */
    std::cout << "Profiles" << std::endl;

    /******** NODE 0 ********/ {
    setting_entry_p_reg node0_set;
    node0_set.i_snum = 0;
    node0_set.Lmin = -300.0;
    node0_set.Lmax = -10.0;
    node0_set.Pmin = 60e5;
    node0_set.Pmax = 80e5;

    std::vector<double> node0_P = {
        70.000000000000000,  69.416666666666671,  68.833333333333343,
        68.250000000000000,  67.666666666666671,  67.083333333333343,
        66.500000000000000,  67.083333333333329,  67.666666666666657,
        68.249999999999986,  68.833333333333329,  69.416666666666657,
        70.000000000000000,  72.333333333333329,  74.666666666666657,
        76.999999999999972,  79.333333333333300,  81.666666666666629,
        84.000000000000000,  81.666666666666657,  79.333333333333329,
        77.000000000000000,  74.666666666666686,  72.333333333333357,
        70.000000000000000 };
        for (auto& p : node0_P)
            p *= 1e5;

    node0_set.Pprofile.reserve(node0_P.size());
    for (int i = 0; i < node0_P.size(); i++)
        node0_set.Pprofile.push_back({double(3600*i), node0_P[i]});
    database::store(st.db, {}, {node0_set});
    }

    /******** NODE 1 ********/ {
    setting_outlet node1_set;
    node1_set.i_snum = 1;
    std::vector<double> node1_L = {
        20.000000000000000,  19.166666666666668,  18.333333333333336,
        17.500000000000004,  16.666666666666671,  15.833333333333337,
        15.000000000000000,  16.666666666666668,  18.333333333333336,
        20.000000000000000,  21.666666666666664,  23.333333333333329,
        25.000000000000000,  23.333333333333336,  21.666666666666671,
        20.000000000000004,  18.333333333333336,  16.666666666666671,
        15.000000000000000,  16.666666666666668,  18.333333333333336,
        20.000000000000000,  21.666666666666664,  23.333333333333329,
        2.5000000e+01 };
    
    node1_set.Lprofile.reserve(node1_L.size());
    for (int i = 0; i < node1_L.size(); i++)
        node1_set.Lprofile.push_back({double(3600*i), node1_L[i]});
    database::store(st.db, {}, {node1_set});
    }

    /******** NODE 5 ********/ {
    setting_exit_l_reg node5_set;
    node5_set.i_snum = 5;
    node5_set.Lmin = -300.0;
    node5_set.Lmax = -10.0;
    node5_set.Pmin = 60e5;
    node5_set.Pmax = 80e5;
    std::vector<double> node5_L = {
        20.000000000000000,  19.333333333333332,  18.666666666666668,
        18.000000000000000,  17.333333333333336,  16.666666666666668,
        16.000000000000000,  16.666666666666668,  17.333333333333336,
        18.000000000000000,  18.666666666666668,  19.333333333333332,
        20.000000000000000,  19.333333333333332,  18.666666666666668,
        18.000000000000000,  17.333333333333336,  16.666666666666668,
        16.000000000000000,  16.666666666666668,  17.333333333333336,
        18.000000000000000,  18.666666666666668,  19.333333333333332,
        2.0000000e+01 };
    
    node5_set.Lprofile.reserve(node5_L.size());
    for (int i = 0; i < node5_L.size(); i++)
        node5_set.Lprofile.push_back({double(3600*i), node5_L[i]});
    database::store(st.db, {}, {node5_set});
    }

    /******** NODE 8 ********/ {
    setting_exit_l_reg node8_set;
    node8_set.i_snum = 8;
    node8_set.Lmin = -300.0;
    node8_set.Lmax = -10.0;
    node8_set.Pmin = 60e5;
    node8_set.Pmax = 80e5;
    std::vector<double> node8_L = {
        50.000000000000000,  52.083333333333336,  54.166666666666671,
        56.250000000000014,  58.333333333333350,  60.416666666666686,
        62.500000000000000,  58.750000000000000,  55.000000000000007,
        51.250000000000007,  47.500000000000007,  43.750000000000014,
        40.000000000000000,  40.000000000000000,  40.000000000000000,
        40.000000000000000,  40.000000000000000,  40.000000000000000,
        40.000000000000000,  41.666666666666671,  43.333333333333336,
        45.000000000000000,  46.666666666666664,  48.333333333333336,
        5.0000000e+01 };
    node8_set.Lprofile.reserve(node8_L.size());
    for (int i = 0; i < node8_L.size(); i++)
        node8_set.Lprofile.push_back({double(3600*i), node8_L[i]});
    database::store(st.db, {}, {node8_set});
    }

    /******** NODE 10 ********/ {
    setting_exit_l_reg node10_set;
    node10_set.i_snum = 10;
    node10_set.Lmin = -300.0;
    node10_set.Lmax = -10.0;
    node10_set.Pmin = 60e5;
    node10_set.Pmax = 80e5;
    std::vector<double> node10_L = {
        15.000000000000000,  14.750000000000000,  14.500000000000002,
        14.250000000000004,  14.000000000000004,  13.750000000000004,
        13.500000000000000,  14.250000000000000,  15.000000000000000,
        15.750000000000000,  16.500000000000000,  17.250000000000004,
        18.000000000000000,  17.499999999999996,  17.000000000000000,
        16.500000000000000,  16.000000000000004,  15.500000000000005,
        15.000000000000000,  15.250000000000000,  15.499999999999998,
        15.749999999999996,  15.999999999999996,  16.249999999999996,
        1.6500000e+01 } ;
    node10_set.Lprofile.reserve(node10_L.size());
    for (int i = 0; i < node10_L.size(); i++)
        node10_set.Lprofile.push_back({double(3600*i), node10_L[i]});
    database::store(st.db, {}, {node10_set});
    }

    /******** NODE 11 ********/ {
    setting_entry_l_reg node11_set;
    node11_set.i_snum = 11;
    node11_set.Lmin = -300.0;
    node11_set.Lmax = -10.0;
    node11_set.Pmin = 60e5;
    node11_set.Pmax = 80e5;
    node11_set.f = 1.0;
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
    node11_set.Pprofile.reserve(node11_P.size());
    node11_set.Lprofile.reserve(node11_L.size());
    assert(node11_set.Lprofile.size() == node11_set.Pprofile.size());
    for (int i = 0; i < node11_L.size(); i++) {
        node11_set.Pprofile.push_back({double(3600*i), node11_P[i]});
        node11_set.Lprofile.push_back({double(3600*i), node11_L[i]});
    }
    database::store(st.db, {}, {node11_set});
    }

    /******** NODE 12 ********/ {
    setting_outlet node12_set;
    node12_set.i_snum = 12;
    std::vector<double> node12_L = {
        10.000000000000000,   9.666666666666666,   9.333333333333334,
         9.000000000000000,   8.666666666666668,   8.333333333333334,
         8.000000000000000,   8.333333333333334,   8.666666666666668,
         9.000000000000000,   9.333333333333334,   9.666666666666666,
        10.000000000000000,   9.833333333333334,   9.666666666666668,
         9.500000000000002,   9.333333333333336,   9.166666666666670,
         9.000000000000000,   9.583333333333334,  10.166666666666666,
        10.750000000000000,  11.333333333333332,  11.916666666666666,
         1.2500000e+01 };
    node12_set.Lprofile.reserve(node12_L.size());
    for (int i = 0; i < node12_L.size(); i++)
        node12_set.Lprofile.push_back({double(3600*i), node12_L[i]});
    database::store(st.db, {}, {node12_set});
    }

    /* Gas mass fractions */
    std::vector<gas_mass_fractions> gmfs;
    for (int i = 0; i < station_types.size(); i++) {
        gas_mass_fractions gmf;
        gmf.i_snum = i;
        for (int ig = 0; ig < NUM_GASES; ig++)
            gmf.fractions[ig] = 0.0;
        gmf.fractions[+gas::CH4] = 1.0;
        gmfs.push_back( std::move(gmf) );
    }
    database::store(st.db, {}, gmfs);

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