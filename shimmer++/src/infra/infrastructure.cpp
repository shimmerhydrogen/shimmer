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

#include <filesystem>
#include "infra/infrastructure.h"
#include "errors.h"

namespace fs = std::filesystem;

namespace shimmer {

static int
populate_entry_p_reg(const std::vector<setting_entry_p_reg>& settings,
    vertex_properties& vp)
{
    auto itor = lookup(settings, vp.i_snum);
    if ( itor == settings.end() ) {
        std::cerr << "Warning: No data for station " << vp.u_snum;
        std::cerr << " (ReMi w/o pressure control). ";
        std::cerr << "Check limits and/or profiles." << std::endl;
        return SHIMMER_MISSING_DATA;
    }
    const setting_entry_p_reg &setting = *itor;
    assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

    auto limits = convert_limits(setting);
    auto Pset = convert_Pprof(setting);
    auto remi = make_station_entry_p_reg(Pset, limits, limits);
    vp.node_station = std::make_unique<multiple_states_station>(remi);          
    return SHIMMER_SUCCESS;
}

static int
populate_entry_l_reg(const std::vector<setting_entry_l_reg>& settings,
    vertex_properties& vp)
{
    auto itor = lookup(settings, vp.i_snum);
    if ( itor == settings.end() ) {
        std::cerr << "Warning: No data for station " << vp.u_snum;
        std::cerr << " (Injection w/ pressure control). ";
        std::cerr << "Check limits and/or profiles."  << std::endl;
        return SHIMMER_MISSING_DATA;
    }
    const setting_entry_l_reg &setting = *itor;
    assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

    auto limits = convert_limits(setting);
    auto Pset = convert_Pprof(setting);
    auto Lset = convert_Lprof(setting);
    auto f = setting.f;
    auto inj_station = make_station_entry_l_reg(f, Pset, Lset,
                                        limits,
                                        limits);
    vp.node_station = std::make_unique<multiple_states_station>(inj_station);
    return SHIMMER_SUCCESS;
}

static int
populate_exit_l_reg(const std::vector<setting_exit_l_reg>& settings,
    vertex_properties& vp)
{
    auto itor = lookup(settings, vp.i_snum);
    if ( itor == settings.end() ) {
        std::cerr << "Warning: No data for station " << vp.u_snum;
        std::cerr << " (consumption w/o pressure control). ";
        std::cerr << "Check limits and/or profiles."  << std::endl;
        return SHIMMER_MISSING_DATA;
    }
    const setting_exit_l_reg &setting = *itor;
    assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

    auto limits = convert_limits(setting);
    auto Lset = convert_Lprof(setting);
    auto consumption = make_station_exit_l_reg(Lset, limits);
    vp.node_station = std::make_unique<one_state_station>(consumption);
    return SHIMMER_SUCCESS;
}

static int
populate_outlet(const std::vector<setting_outlet>& settings,
    vertex_properties& vp)
{
    auto itor = lookup(settings, vp.i_snum);
    if ( itor == settings.end() ) {
        std::cout << "Warning: No data for station " << vp.u_snum;
        std::cout << " (Outlet)" << std::endl;
        return SHIMMER_MISSING_DATA;
    }
    const setting_outlet &setting = *itor;
    assert( setting.i_snum == vp.i_snum );

    auto Lset = convert_Lprof(setting);
    auto exit_station = priv::make_station_outlet(Lset);
    vp.node_station = std::make_unique<one_state_station>(exit_station);
    return SHIMMER_SUCCESS;
}

static int
populate_type_dependent_station_data(const infrastructure& infra,
    vertex_properties& vp)
{
    int err = SHIMMER_SUCCESS;
    switch (vp.type) {
        case station_type::JUNCTION:
            vp.node_station = std::make_unique<junction>();
            break;

        case station_type::FICTITIOUS_JUNCTION:
            vp.node_station = std::make_unique<fictitious_junction>();
            break;

        case(station_type::ENTRY_P_REG):
            err = populate_entry_p_reg(infra.settings_entry_p_reg, vp);
            break;
            
        case(station_type::ENTRY_L_REG):
            err = populate_entry_l_reg(infra.settings_entry_l_reg, vp);
            break;

        case(station_type::EXIT_L_REG):
            err = populate_exit_l_reg(infra.settings_exit_l_reg, vp);
            break;
        
        
        case(station_type::PRIVATE_OUTLET):
            err = populate_outlet(infra.settings_outlet, vp);
            break;
            
        default:
            std::cerr << "WARNING: Unhandled station type " << +vp.type;
            std::cerr << " while loading type-dependent station data" << std::endl;
            return SHIMMER_INVALID_DATA;
    }

    return err;  
}

static int
load_stations(sqlite3 *db, infrastructure& infra)
{              
    int rc;
    std::string q = "SELECT * FROM stations";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    infra.s_u2vd.resize( infra.s_u2i.size() );
    infra.s_i2vd.resize( infra.s_i2u.size() );

    int incomplete_stations = 0;

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        vertex_properties vp;
        
        /* Station number, convert from user to internal */
        int s_un = sqlite3_column_int(stmt, 0);
        int s_in = infra.s_u2i.at(s_un).value();
        vp.u_snum = s_un;
        vp.i_snum = s_in;
        /* Station name and type */
        vp.name = (char *) sqlite3_column_text(stmt, 1);
        vp.type = static_cast<station_type>(sqlite3_column_int(stmt,2));

        if ( s_in >= infra.molar_fractions.size() ) {
            std::cerr << "No mass fractions for station " << s_un << std::endl;
            return SHIMMER_MISSING_DATA;
        }
        const auto& mfs = infra.molar_fractions.at(s_in).fractions;
        vp.gas_mixture = vector_t::Zero(mfs.size());
        assert(vp.gas_mixture.size() == NUM_GASES);
        std::copy(mfs.begin(), mfs.end(), vp.gas_mixture.begin());

        /* Station location */
        vp.height = sqlite3_column_double(stmt, 3);
        vp.latitude = sqlite3_column_double(stmt, 4);
        vp.longitude = sqlite3_column_double(stmt, 5);
        
        if (populate_type_dependent_station_data(infra, vp) != SHIMMER_SUCCESS) {
            incomplete_stations++;
        }
        
        vertex_descriptor vtxd = boost::add_vertex(infra.graph);
        /* keep track of the mapping from user/internal numbering to boost
         * vertex descriptors. this should be done with boost::labeled_graph,
         * however it does not like unique pointers, so here we are. */
        infra.s_u2vd.at(s_un) = vtxd;
        infra.s_i2vd.at(s_in) = vtxd;
        
        /* Move the newly created vertex_properties to the graph */
        infra.graph[vtxd] = std::move(vp);
    }

    sqlite3_finalize(stmt);

    if (incomplete_stations > 0)
        return SHIMMER_MISSING_DATA;

    return SHIMMER_SUCCESS;
}

static int
store_stations(sqlite3 *db, infrastructure& infra)
{
    int rc;
    std::string q = "INSERT INTO stations VALUES "
        "(?, ?, ?, ?, ?, ?)";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);
    auto [vbegin, vend] = boost::vertices(infra.graph);
    for (auto vitor = vbegin; vitor != vend; vitor++)
    {
        auto& vp = infra.graph[*vitor];
        assert( vp.u_snum == convert_i2u(infra.s_i2u, vp.i_snum) );
        rc = sqlite3_bind_int(stmt, 1, vp.u_snum);
        rc = sqlite3_bind_text(stmt, 2, vp.name.c_str(), vp.name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 3, +vp.type);
        rc = sqlite3_bind_double(stmt, 4, vp.height);
        rc = sqlite3_bind_double(stmt, 5, vp.latitude);
        rc = sqlite3_bind_double(stmt, 6, vp.longitude);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    return SHIMMER_SUCCESS;
}

static int
populate_type_dependent_pipe_data(infrastructure& infra,
    edge_properties& ep, int i_from, int i_to)
{
    using edge_constraint_t = edge_station::control::constraint_type;
    using external_t = edge_station::external_type; 
    using user_limits_t = std::unordered_map<external_t, std::pair<edge_constraint_t, double>>;

    switch (ep.type) {
        case pipe_type::PIPE: {
            auto sitor = lookup(infra.settings_pipe, i_from, i_to);
            if (sitor == infra.settings_pipe.end()) {
                auto u_from = infra.s_i2u.at(i_from);
                auto u_to = infra.s_i2u.at(i_to);
                std::cerr << "WARNING: no data for pipe (" << u_from << ", ";
                std::cerr << u_to << ")" << std::endl;
                return SHIMMER_MISSING_DATA;
            }

            ep.friction_factor = (*sitor).roughness;
            ep.diameter = (*sitor).diameter;
            ep.length = (*sitor).length;
            break;
        }
        case pipe_type::COMPR_STAT: {
            auto sitor = lookup(infra.settings_compr_stat, i_from, i_to);
            if (sitor == infra.settings_compr_stat.end()) {
                auto u_from = infra.s_i2u.at(i_from);
                auto u_to = infra.s_i2u.at(i_to);
                std::cerr << "WARNING: no data for compressor (";
                std::cerr << u_from << ", " << u_to << ")" << std::endl;
                return SHIMMER_MISSING_DATA;
            }
            auto& setting = *sitor;

            user_limits_t user_limits;
            
            user_limits[external_t::P_OUT_MAX] =
                std::pair(edge_constraint_t::LOWER_EQUAL, setting.max_outpress);
            user_limits[external_t::P_IN_MIN] =
                std::pair(edge_constraint_t::GREATER_EQUAL, setting.min_inpress);
            user_limits[external_t::BETA_MAX] =
                std::make_pair(edge_constraint_t::LOWER_EQUAL, setting.max_ratio);
            user_limits[external_t::BETA_MIN] =
                std::make_pair(edge_constraint_t::GREATER_EQUAL, setting.min_ratio);
            user_limits[external_t::FLUX_MAX] =
                std::make_pair(edge_constraint_t::LOWER_EQUAL, setting.max_massflow);
            user_limits[external_t::PWD_NOMINAL] =
                std::make_pair(edge_constraint_t::LOWER_EQUAL, setting.max_power);

            std::vector<std::pair<compressor_mode,double>> mode_type_vec;
            mode_type_vec.resize( setting.profile.size() );

            /* settings.profile.size() must have as many entries as required timesteps */

            std::transform(setting.profile.begin(), setting.profile.end(),
                mode_type_vec.begin(),
                [](const compressor_profile_sample& cps){ return cps.value_bymode(); }
            );

            auto num_steps = setting.profile.size();
            std::vector<bool> activate_history (num_steps, true); 

            auto comp = edge_station::make_compressor(num_steps, 
                                                      setting.ramp_coeff,
                                                      setting.efficiency, 
                                                      activate_history,
                                                      mode_type_vec,
                                                      user_limits);

            ep.pipe_station = std::make_shared<edge_station::compressor>(comp);

            break;
        }
        case pipe_type::RED_STAT: {
            auto sitor = lookup(infra.settings_red_stat, i_from, i_to);
            if (sitor == infra.settings_red_stat.end()) {
                auto u_from = infra.s_i2u.at(i_from);
                auto u_to = infra.s_i2u.at(i_to);
                std::cerr << "WARNING: no data for reduction station (";
                std::cerr << u_from << ", " << u_to << ")" << std::endl;
                return SHIMMER_MISSING_DATA;
            }

            auto& setting = *sitor;
            /* extract & populate from settin */

            std::cerr << "RED_STAT not implemented" << std::endl;
            return SHIMMER_INVALID_DATA;
            break;
        }

        case pipe_type::VALVE: {
            std::cerr << "VALVE not implemented" << std::endl;
            return SHIMMER_INVALID_DATA;
            break;
        }

        default:
            std::cerr << "WARNING: invalid pipe type " << +ep.type << std::endl;
            return SHIMMER_INVALID_DATA;
    }

    return SHIMMER_SUCCESS;
}

static int
load_pipelines(sqlite3 *db, infrastructure& infra)
{
    int rc;
    std::string q = "SELECT * FROM pipelines";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    int branch_num = 0;
    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int from_un = sqlite3_column_int(stmt, 1);
        int to_un = sqlite3_column_int(stmt, 2);
        auto from_vtx = infra.s_u2vd.at(from_un).value();
        auto to_vtx = infra.s_u2vd.at(to_un).value();

        auto from_in = infra.s_u2i.at(from_un).value();
        auto to_in = infra.s_u2i.at(to_un).value();
        
        /* XXX: NAME CURRENTLY NOT USED - CANNOT HAVE
         *      MULTIPLE PIPES BETWEEN THE SAME TWO NODES
         */
        int type = sqlite3_column_int(stmt, 3);
        
        edge_properties ep;
        ep.name = name;
        ep.i_sfrom = from_in;
        ep.i_sto = to_in;
        ep.type = static_cast<pipe_type>(type);
        ep.branch_num = branch_num++;

        populate_type_dependent_pipe_data(infra, ep, from_in, to_in);

        boost::add_edge(from_vtx, to_vtx, ep, infra.graph);
    }
    //npipes_ = branch_num;

    sqlite3_finalize(stmt);

    return SHIMMER_SUCCESS;
}

static int
store_pipelines(sqlite3 *db, infrastructure& infra)
{
    int rc;
    std::string q = "INSERT INTO pipelines VALUES "
        "(?, ?, ?, ?)";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);
    auto [ebegin, eend] = boost::edges(infra.graph);
    for (auto eitor = ebegin; eitor != eend; eitor++)
    {
        auto& ep = infra.graph[*eitor];
        auto u_from = convert_i2u(infra.s_i2u, ep.i_sfrom);
        auto u_to = convert_i2u(infra.s_i2u, ep.i_sto);
        rc = sqlite3_bind_text(stmt, 1, ep.name.c_str(), ep.name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, u_from);
        rc = sqlite3_bind_int(stmt, 3, u_to);
        rc = sqlite3_bind_int(stmt, 4, +ep.type);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    return SHIMMER_SUCCESS;
}

static int
load_station_settings(sqlite3 *db, infrastructure& infra)
{
    int errors = 0;
    if (SHIMMER_SUCCESS != database::renumber_stations(db, infra.s_u2i, infra.s_i2u))
        errors++;
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_entry_p_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_entry_l_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_exit_l_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_outlet))
        errors++;

    if (errors > 0)
        return SHIMMER_INVALID_DATA;

    std::cout << "-- Non-junction station summary --" << std::endl;
    std::cout << "  Type Entry P:" << std::endl;
    for (auto& s : infra.settings_entry_p_reg)
        std::cout << "    " << s << std::endl;
    std::cout << "  Type Entry L:" << std::endl;
    for (auto& s : infra.settings_entry_l_reg)
        std::cout << "    " << s << std::endl;
    std::cout << "  Type Exit L:" << std::endl;
    for (auto& s : infra.settings_exit_l_reg)
        std::cout << "    " << s << std::endl;
    std::cout << "  Type outlet:" << std::endl;
    for (auto& s : infra.settings_outlet)
        std::cout << "    " << s << std::endl;

    return SHIMMER_SUCCESS;
}

static int
store_station_settings(sqlite3 *db, infrastructure& infra)
{
    int errors = 0;
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_entry_p_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_entry_l_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_exit_l_reg))
        errors++;
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_outlet))
        errors++;

    if (errors > 0)
        return SHIMMER_INVALID_DATA;

    return SHIMMER_SUCCESS;
}

static int
load_gas_molar_fractions(sqlite3 *db, infrastructure& infra)
{
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.molar_fractions)) {
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    if (infra.molar_fractions.size() != infra.s_i2u.size()) {
        std::cerr << "Incorrect number of mass fractions in the DB. ";
        std::cerr << "There should be as many mass fractions as stations.";
        std::cerr << std::endl;
        return SHIMMER_INVALID_DATA;
    }
    
    for (int i = 0; i < infra.s_i2u.size(); i++) {
        if (infra.molar_fractions[i].i_snum != i) {
            std::cerr << "Incoherent mass fraction data. Check that ";
            std::cerr << "there is 1:1 correspondence between stations ";
            std::cerr << "and mass fraction data" << std::endl;
            return SHIMMER_INVALID_DATA;
        }
    }
    
    return SHIMMER_SUCCESS;
}

static int
count_edges_bytype(const infrastructure& infra, pipe_type type)
{
    auto edge_range = boost::edges(infra.graph);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    size_t found_pipes = 0;
    for (auto itor = begin; itor != end; itor++) {
        if (infra.graph[*itor].type == type)
            found_pipes++;
    }
    return found_pipes;
}

static int
check_pipeline_data_consistency(infrastructure& infra)
{
    /* pipes */ {
    int num_edges = count_edges_bytype(infra, pipe_type::PIPE);
    int num_settings = infra.settings_pipe.size();

    if (num_edges != num_settings) {
        std::cerr << "Error: In the database there are " << num_edges;
        std::cerr << " pipes but " << num_settings << " settings.";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    }

    /* compressors */ {
    int num_edges = count_edges_bytype(infra, pipe_type::COMPR_STAT);
    int num_settings = infra.settings_compr_stat.size();

    if (num_edges != num_settings) {
        std::cerr << "Error: In the database there are " << num_edges;
        std::cerr << " compressors but " << num_settings << " settings.";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    }

    return SHIMMER_SUCCESS;
}

int load(const std::string& db_filename, infrastructure& infra)
{
    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != load_station_settings(db, infra)) {
        std::cerr << "Problems detected while loading station settings. ";
        std::cerr << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != load_gas_molar_fractions(db, infra)) {
        std::cerr << "Problems detected while loading gas molar fractions. ";
        std::cerr << std::endl;
        goto load_fail;
    }
    
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_pipe)) {
        std::cerr << "Problems detected while loading pipe settings";
        std::cerr << std::endl;
        goto load_fail;
    };
    
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_compr_stat)) {
        std::cerr << "Problems detected while loading compressor settings";
        std::cerr << std::endl;
        goto load_fail;
    };

    if (SHIMMER_SUCCESS != load_stations(db, infra)) {
        std::cerr << "Problem detected while loading stations";
        std::cerr << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != load_pipelines(db, infra)) {
        std::cerr << "Problem detected while loading pipelines";
        std::cerr << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != check_pipeline_data_consistency(infra)) {
        std::cerr << "Inconsistencies detected in pipeline data";
        std::cerr << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.sics) ) {
        std::cerr << "Problem detected while loading initial condition for stations";
        std::cerr << std::endl;
        goto load_fail;
    }

    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.pics) ) {
        std::cerr << "Problem detected while loading initial condition for pipelines";
        std::cerr << std::endl;
        goto load_fail;
    }

    sqlite3_close(db);
    return SHIMMER_SUCCESS;

load_fail:
    sqlite3_close(db);
    return SHIMMER_DATABASE_PROBLEM;
}

int store(const std::string& db_filename, infrastructure& infra)
{
    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != store_station_settings(db, infra)) {
        std::cerr << "Problems detected while storing station settings. ";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    //if (SHIMMER_SUCCESS != load_gas_molar_fractions(db, infra)) {
    //    std::cerr << "Problems detected while loading gas mass fractions. ";
    //    std::cerr << std::endl;
    //    return SHIMMER_DATABASE_PROBLEM;
    //}

    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_pipe)) {
        std::cerr << "Problems detected while storing pipe settings";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    };
    
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.settings_compr_stat)) {
        std::cerr << "Problems detected while storing compressor settings";
        std::cerr << std::endl;
        //return SHIMMER_DATABASE_PROBLEM;
    };
    
    if (SHIMMER_SUCCESS != store_stations(db, infra)) {
        std::cerr << "Problem detected while storing stations";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != store_pipelines(db, infra)) {
        std::cerr << "Problem detected while storing pipelines";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.sics) ) {
        std::cerr << "Problem detected while storing initial condition for stations";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    if (SHIMMER_SUCCESS != database::store(db, infra.s_i2u, infra.pics) ) {
        std::cerr << "Problem detected while storing initial condition for pipelines";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}


int num_stations(const infrastructure& infra)
{
    assert(infra.s_i2u.size() == num_vertices(infra.graph));
    return infra.s_i2u.size();
}

int num_pipes(const infrastructure& infra)
{
    return num_edges(infra.graph);
}

variable
initial_guess(const infrastructure& infra)
{
    int nstations = num_stations(infra);
    int npipes = num_pipes(infra);

    vector_t P = vector_t::Zero(nstations);
    vector_t L = vector_t::Zero(nstations);
    vector_t G = vector_t::Zero(npipes);

    for (const auto& sic : infra.sics) {
        if (sic.i_snum >= nstations) {
            throw std::logic_error("station index out of bounds");
        }

        if (sic.init_P == 0 or sic.init_L == 0) {
            std::cout << "Warning: zero initial guess for station ";
            std::cout << infra.s_i2u.at(sic.i_snum) << std::endl;
        }

        P(sic.i_snum) = sic.init_P;
        L(sic.i_snum) = sic.init_L;
    }

    for (size_t i = 0; i < infra.pics.size(); i++) {
        assert(i < npipes);
        const auto& pic = infra.pics[i];

        if (pic.init_G == 0) {
            std::cout << "Warning: zero initial guess for pipe ";
            std::cout << infra.s_i2u.at(pic.i_sfrom) << "-";
            std::cout << infra.s_i2u.at(pic.i_sto) << std::endl;
        }

        G(i) = pic.init_G;
    }

    return variable(P, G, L);
}

int save_pressures(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& pressures)
{
    assert(pressures.cols() == infra.s_i2u.size());

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_station_pressures VALUES (?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_station_pressures", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < pressures.rows(); ts++) {
        for (int i_snum = 0; i_snum < pressures.cols(); i_snum++) {
            rc = sqlite3_bind_int(stmt, 1, convert_i2u(infra.s_i2u, i_snum));
            rc = sqlite3_bind_int(stmt, 2, ts);
            rc = sqlite3_bind_double(stmt, 3, pressures(ts,i_snum));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

int save_flowrates(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& flowrates)
{
    assert(flowrates.cols() == num_pipes(infra));

    auto edge_range = boost::edges(infra.graph);
    auto edge_begin = edge_range.first;
    auto edge_end = edge_range.second;

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_pipe_flowrates VALUES (?, ?, ?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_pipe_flowrates", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < flowrates.rows(); ts++) {
        for (int branch_num = 0; branch_num < flowrates.cols(); branch_num++) {
            auto edge_itor = std::next(edge_begin, branch_num);
            auto ep = infra.graph[*edge_itor];
            rc = sqlite3_bind_text(stmt, 1, ep.name.c_str(), -1, nullptr);
            rc = sqlite3_bind_int(stmt, 2, convert_i2u(infra.s_i2u, ep.i_sfrom));
            rc = sqlite3_bind_int(stmt, 3, convert_i2u(infra.s_i2u, ep.i_sto));
            rc = sqlite3_bind_int(stmt, 4, ts);
            rc = sqlite3_bind_double(stmt, 5, flowrates(ts, branch_num));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

int save_flowrates_stations(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& flowrates_stations)
{
    assert(flowrates_stations.cols() == infra.s_i2u.size());

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_station_flowrates VALUES (?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_station_flowrates", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < flowrates_stations.rows(); ts++) {
        for (int i_snum = 0; i_snum < flowrates_stations.cols(); i_snum++) {
            rc = sqlite3_bind_int(stmt, 1, convert_i2u(infra.s_i2u, i_snum));
            rc = sqlite3_bind_int(stmt, 2, ts);
            rc = sqlite3_bind_double(stmt, 3, flowrates_stations(ts,i_snum));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

int save_velocities(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& vels)
{
    assert(vels.cols() == num_pipes(infra));

    auto edge_range = boost::edges(infra.graph);
    auto edge_begin = edge_range.first;
    auto edge_end = edge_range.second;

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_pipe_velocities VALUES (?, ?, ?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_pipe_velocities", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (int ts = 0; ts < vels.rows(); ts++) {
        std::cout << ts << ": ";
        for (int branch_num = 0; branch_num < vels.cols(); branch_num++) {
            std::cout << branch_num << " ";
            auto edge_itor = std::next(edge_begin, branch_num);
            auto ep = infra.graph[*edge_itor];
            rc = sqlite3_bind_text(stmt, 1, ep.name.c_str(), -1, nullptr);
            rc = sqlite3_bind_int(stmt, 2, convert_i2u(infra.s_i2u, ep.i_sfrom));
            rc = sqlite3_bind_int(stmt, 3, convert_i2u(infra.s_i2u, ep.i_sto));
            rc = sqlite3_bind_int(stmt, 4, ts);
            rc = sqlite3_bind_double(stmt, 5, vels(ts, branch_num));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
        std::cout << std::endl;
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

int save_molar_fractions(const std::string& db_filename, const infrastructure& infra, std::vector<matrix_t>& molar_fractions)
{
    assert(boost::num_vertices(infra.graph) == infra.s_i2u.size());

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_station_molarfrac VALUES (?, ?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_station_molarfrac", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < molar_fractions.size(); ts++) {
        for (int i_snum = 0; i_snum < molar_fractions[ts].rows(); i_snum++) {
            for (int comp = 0; comp < molar_fractions[ts].cols(); comp++) {
                rc = sqlite3_bind_int(stmt, 1, convert_i2u(infra.s_i2u, i_snum));
                rc = sqlite3_bind_int(stmt, 2, ts);
                rc = sqlite3_bind_int(stmt, 3, comp);
                rc = sqlite3_bind_double(stmt, 4, molar_fractions[ts](i_snum, comp));
                rc = sqlite3_step(stmt);
                rc = sqlite3_clear_bindings(stmt);
                rc = sqlite3_reset(stmt);
            }
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}


static void
transfer_original_stations(const infrastructure& infrain,
    infrastructure& infraout)
{
    infraout.s_u2i = infrain.s_u2i;
    infraout.s_i2u = infrain.s_i2u;
    infraout.s_u2vd.resize( infrain.s_u2vd.size() );
    infraout.s_i2vd.resize( infrain.s_i2vd.size() );

    int incomplete_stations = 0;
    auto [vbegin, vend] = vertices(infrain.graph);
    for (auto vitor = vbegin; vitor != vend; vitor++) {
        const vertex_properties& invp = infrain.graph[*vitor];
        vertex_properties outvp;
        outvp.u_snum = invp.u_snum;
        outvp.i_snum = invp.i_snum;
        outvp.name = invp.name;
        outvp.type = invp.type;
        outvp.gas_mixture = invp.gas_mixture;
        outvp.height = invp.height;
        outvp.latitude = invp.latitude;
        outvp.longitude = invp.longitude;

        if (populate_type_dependent_station_data(infraout, outvp) != SHIMMER_SUCCESS) {
            incomplete_stations++;
        }

        if (incomplete_stations) {
            std::cout << "transfer_original_stations(): incomplete data" << std::endl;
        }
        
        vertex_descriptor vtxd = boost::add_vertex(infraout.graph);
        infraout.s_u2vd.at(outvp.u_snum) = vtxd;
        infraout.s_i2vd.at(outvp.i_snum) = vtxd;
        
        infraout.graph[vtxd] = std::move(outvp);
    }
}

static void
discretize_pipes(const infrastructure& infrain,
    infrastructure& infraout, double dx)
{
    int unum_max = 0;
    auto [vbegin, vend] = boost::vertices(infrain.graph);
    for (auto vtxitor = vbegin; vtxitor != vend; vtxitor++) {
        auto& vp = infrain.graph[*vtxitor];
        unum_max = std::max(vp.u_snum, unum_max);
    }
    assert(unum_max+1 == infraout.s_u2i.size());

    size_t fict_station_counter = 0;
    size_t fict_station_ubase = unum_max+1;
    size_t fict_station_ibase = num_stations(infrain);

    int branch_num = 0;
    auto [ebegin, eend] = boost::edges(infrain.graph);

    // Pipe elements
    for (auto edgeitor = ebegin; edgeitor != eend; edgeitor++)
    {
        auto pipe = infrain.graph[*edgeitor];
        if (pipe.type != pipe_type::PIPE) {
            continue;
        }

        auto s = boost::source(*edgeitor, infrain.graph);
        auto t = boost::target(*edgeitor, infrain.graph);

        auto i_from = infrain.graph[s].i_snum;
        auto u_from = infrain.graph[s].u_snum;
        auto lat_from = infrain.graph[s].latitude;
        auto lon_from = infrain.graph[s].longitude;

        auto i_to = infrain.graph[t].i_snum;
        auto u_to = infrain.graph[t].u_snum;
        auto lat_to = infrain.graph[t].latitude;
        auto lon_to = infrain.graph[t].longitude;


        auto settingitor = lookup(infrain.settings_pipe, i_from, i_to);
        if (settingitor == infrain.settings_pipe.end()) {
            throw std::logic_error("pipe not found");
        }
        auto in_setting = *settingitor;

        double numfrags = 1.0;
        if (in_setting.ref_nsegs == 0) {
            numfrags = std::ceil(std::abs(in_setting.length/dx));
            std::cout << "Pipe " << u_from << " -> " << u_to << ": split in ";
            std::cout << numfrags << " segments using global dx" << std::endl;
        } else {
            numfrags = in_setting.ref_nsegs;
            std::cout << "Pipe " << u_from << " -> " << u_to << ": split in ";
            std::cout << numfrags << " segments using database setting" << std::endl;
        }
        
        auto fraglen = in_setting.length/numfrags;
        auto dlat = (lat_to - lat_from)/numfrags;
        auto dlon = (lon_to - lon_from)/numfrags;

        std::string nname = "fict_(" + std::to_string(u_from) + ","
            + std::to_string(u_to) + ")_";

        assert(numfrags > 0);
        std::vector<int> discrnodes(numfrags+1);

        auto from_ic = lookup(infrain.sics, i_from);
        auto to_ic = lookup(infrain.sics, i_to);
        auto pipe_ic = lookup(infrain.pics, i_from, i_to);

        auto interp = [](double x, double x2, double q1, double q2) {
            return x*(q2 - q1)/x2 + q1;
        };

        // Add discretization nodes for each original pipe
        discrnodes[0] = i_from;
        vector_t methane = vector_t::Zero(NUM_GASES); 
        methane(GAS_TYPE::CH4) = 1.0;
        //std:: cout << from_ic->init_P << " -> " << from_ic->init_L << std::endl;
        for (int i = 1; i < numfrags; i++) {
            auto inum = fict_station_ibase + fict_station_counter;
            auto unum = fict_station_ubase + fict_station_counter;
            discrnodes[i] = inum;
            vertex_properties outvp;
            outvp.u_snum = unum;
            outvp.i_snum = inum;
            outvp.name = nname + std::to_string(unum);
            outvp.type = station_type::FICTITIOUS_JUNCTION;
            outvp.gas_mixture = methane;

            outvp.latitude = lat_from + i*dlat;
            outvp.longitude = lon_from + i*dlon;

            vertex_descriptor vtxd = boost::add_vertex(infraout.graph);
            assert(infraout.s_i2vd.size() == outvp.i_snum);
            infraout.s_i2vd.push_back(vtxd);
            populate_type_dependent_station_data(infraout, outvp);
            infraout.graph[vtxd] = std::move(outvp);
            infraout.s_i2u.push_back(unum);
            
            auto x = i*fraglen;
            station_initial_condition sic;
            sic.i_snum = inum;
            sic.init_P = interp(x, in_setting.length, from_ic->init_P, to_ic->init_P);
            sic.init_L = 0.0;//interp(x, in_setting.length, from_ic->init_L, to_ic->init_L);
            //std::cout << x << " -> " << sic.init_P << " -> " << sic.init_L << std::endl;
            infraout.sics.push_back(sic);

            fict_station_counter++;
        }
        discrnodes[numfrags] = i_to;
        //std:: cout << to_ic->init_P << " -> " << to_ic->init_L << std::endl;

        // Add discretization pipes for each original pipe
        std::vector<int> discrpipes(numfrags);
        for (int i = 1; i < discrnodes.size(); i++) {
            setting_pipe out_setting;
            out_setting.name = pipe.name + "_seg_" + std::to_string(i);
            out_setting.i_sfrom = discrnodes[i-1];
            out_setting.i_sto = discrnodes[i];
            out_setting.length = fraglen;
            out_setting.diameter = in_setting.diameter;
            out_setting.roughness = in_setting.roughness;
            out_setting.ref_nsegs = 0;
            infraout.settings_pipe.push_back(out_setting);
            
            edge_properties newnp;
            newnp.type = pipe_type::PIPE;
            newnp.branch_num = branch_num;
            newnp.length = pipe.length;
            newnp.diameter = pipe.diameter;
            newnp.friction_factor = pipe.friction_factor;
            newnp.name = out_setting.name;
            newnp.i_sfrom = out_setting.i_sfrom;
            newnp.i_sto = out_setting.i_sto;
            auto from_vtx = infraout.s_i2vd[newnp.i_sfrom];
            auto to_vtx = infraout.s_i2vd[newnp.i_sto];           
            auto [ed, is_added] = boost::add_edge(from_vtx, to_vtx, newnp, infraout.graph);
            infraout.p_i2ed.push_back(ed);

            pipe_initial_condition pic;
            pic.i_sfrom = discrnodes[i-1];
            pic.i_sto = discrnodes[i];
            pic.init_G = pipe_ic->init_G;
            infraout.pics.push_back(pic);

            discrpipes.at(i-1) = branch_num;

            branch_num++;
        }

        pipe_discretization discr;
        discr.parent_ifrom = i_from;
        discr.parent_ito = i_to;
        discr.dx = fraglen;
        discr.nodelist = std::move(discrnodes);
        discr.pipelist = std::move(discrpipes); 

        infraout.pipe_discretizations.push_back( std::move(discr) );
    }

    int unum_max_2 = infraout.s_i2u.back();
    infraout.s_u2i.resize( unum_max_2+1 );
    infraout.s_u2vd.resize( unum_max_2+1 );
    for (size_t i = 0; i < infraout.s_i2u.size(); i++) {
        infraout.s_u2i[ infraout.s_i2u[i] ] = infraout.s_i2u[i];
        infraout.s_u2vd[ infraout.s_i2u[i] ] = infraout.s_i2vd[i];

    }

    // Non-pipe elements
    for (auto edgeitor = ebegin; edgeitor != eend; edgeitor++)
    {
        auto pipe = infrain.graph[*edgeitor];
        if (pipe.type == pipe_type::PIPE) {
            continue;
        }

        edge_properties newnp;
        newnp.type = pipe.type;
        newnp.branch_num = branch_num;
        newnp.length = pipe.length;
        newnp.diameter = pipe.diameter;
        newnp.friction_factor = pipe.friction_factor;
        newnp.name = pipe.name;
        newnp.i_sfrom = pipe.i_sfrom;
        newnp.i_sto = pipe.i_sto;
        populate_type_dependent_pipe_data(infraout, newnp, newnp.i_sfrom, newnp.i_sto);
        auto from_vtx = infraout.s_i2vd[newnp.i_sfrom];
        auto to_vtx = infraout.s_i2vd[newnp.i_sto];
        //std::cout << from_vtx << " " << to_vtx << std::endl;       
        auto [ed, is_added] = boost::add_edge(from_vtx, to_vtx, newnp, infraout.graph);
        infraout.p_i2ed.push_back(ed);

        branch_num++;
    }

    std::sort(infraout.pipe_discretizations.begin(),
        infraout.pipe_discretizations.end());
    std::sort(infraout.sics.begin(), infraout.sics.end());
    std::sort(infraout.pics.begin(), infraout.pics.end());
}

int
refine_pipes(const infrastructure& infrain,
    infrastructure& infraout, double dx)
{
    /* load_station_settings() */
    infraout.settings_outlet = infrain.settings_outlet;
    infraout.settings_entry_p_reg = infrain.settings_entry_p_reg;
    infraout.settings_entry_l_reg = infrain.settings_entry_l_reg;
    infraout.settings_exit_l_reg = infrain.settings_exit_l_reg;
    infraout.sics = infrain.sics;
    infraout.num_original_stations = num_stations(infrain);
    infraout.num_original_pipes = num_pipes(infrain);

    //infraout.pics = infrain.pics;
    /* This has to be recomputed:
     * load_gas_molar_fractions() */

    /* load_stations() */
    transfer_original_stations(infrain, infraout);

    /* This has to be recomputed: 
     * load(db, infra.s_u2i, infra.settings_pipe) */
    
    discretize_pipes(infrain, infraout, dx);

    assert(infraout.sics.size() == num_vertices(infraout.graph));
    std::cout << infraout.pics.size() << " -- " << boost::num_edges(infraout.graph) << std::endl;
    assert(infraout.pics.size() == num_edges(infraout.graph));

    /* This remains the same:
     * load(db, infra.s_u2i, infra.settings_compr_stat) */
    infraout.settings_compr_stat = infrain.settings_compr_stat;
    infraout.settings_red_stat = infrain.settings_red_stat;

    

    /* load_pipelines() */
    


    

    return 0;
}

config::config()
{
    database = "";
    steps = 1;
    dt_std = 1;
    dt = 3600;
    temperature = 293.15;
    tol_std = 1e-14;
    tol = 1e-4;
}

int initdb(const std::string& db_filename)
{
    static const size_t BUFSIZE = 1024;

    char buf[BUFSIZE];
    FILE *fh;

    if ( (fh = fopen("../../sqlite/shimmer.sql", "r")) )
        goto foundok;

    if ( (fh = fopen("shimmer.sql", "r")) )
        goto foundok;

    std::cerr << "Can't find database schema file. ";
    std::cerr << "Can't proceed." << std::endl;
    return -1;

foundok:
    unlink(db_filename.c_str());
    std::string popenstr = "sqlite3 " + db_filename;

    FILE *ph;
    ph = popen(popenstr.c_str(), "w");
    if (not ph) {
        std::cerr << "Failed to populate DB: system error" << std::endl;
        return -1;
    }

    while (!feof(fh)) {
        size_t nread = fread(buf, 1, BUFSIZE, fh);
        fwrite(buf, 1, nread, ph);
    }

    pclose(fh);
    int sqlite_ret = pclose(ph);
    
    return sqlite_ret;
}

} // namespace shimmer
