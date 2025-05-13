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

#include "infra/infrastructure.h"
#include "errors.h"

namespace shimmer {

static int
populate_entry_p_reg(const std::vector<setting_entry_p_reg>& settings,
    vertex_properties& vp)
{
    auto itor = lookup(settings, vp.i_snum);
    if ( itor == settings.end() ) {
        std::cout << "Warning: No data for station " << vp.u_snum;
        std::cout << " (ReMi w/o pressure control)" << std::endl;
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
        std::cout << "Warning: No data for station " << vp.u_snum;
        std::cout << " (Injection w/ pressure control)" << std::endl;
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
        std::cout << "Warning: No data for station " << vp.u_snum;
        std::cout << " (consumption w/o pressure control)" << std::endl;
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

        if ( s_in >= infra.mass_fractions.size() ) {
            std::cerr << "No mass fractions for station " << s_un << std::endl;
            return SHIMMER_MISSING_DATA;
        }
        const auto& mfs = infra.mass_fractions.at(s_in).fractions;
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

            std::transform(setting.profile.begin(), setting.profile.end(),
                mode_type_vec.begin(),
                [](const compressor_profile_sample& cps){ return cps.value_bymode(); }
            );

            auto num_steps = 0;
            std::vector<bool> activate_history ( num_steps, true); 

            auto comp = edge_station::make_compressor(setting.ramp_coeff,
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
load_gas_mass_fractions(sqlite3 *db, infrastructure& infra)
{
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.mass_fractions)) {
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    if (infra.mass_fractions.size() != infra.s_i2u.size()) {
        std::cerr << "Incorrect number of mass fractions in the DB. ";
        std::cerr << "There should be as many mass fractions as stations.";
        std::cerr << std::endl;
        return SHIMMER_INVALID_DATA;
    }
    
    for (int i = 0; i < infra.s_i2u.size(); i++) {
        if (infra.mass_fractions[i].i_snum != i) {
            std::cerr << "Incoherent mass fraction data. Check that ";
            std::cerr << "there is 1:1 correspondence between stations ";
            std::cerr << "and mass fraction data" << std::endl;
            return SHIMMER_INVALID_DATA;
        }
    }
    
    return SHIMMER_SUCCESS;
}

int load(const std::string db_filename, infrastructure& infra)
{
    sqlite3 *db;
    int rc = sqlite3_open(db_filename.c_str(), &db);
    if(rc) {
        std::cerr << "Can't open database: ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != load_station_settings(db, infra)) {
        std::cerr << "Problems detected while loading station settings. ";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != load_gas_mass_fractions(db, infra)) {
        std::cerr << "Problems detected while loading gas mass fractions. ";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_pipe)) {
        std::cerr << "Problems detected while loading pipe settings";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    };

    if (SHIMMER_SUCCESS != load_stations(db, infra)) {
        std::cerr << "Problem detected while loading stations";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != load_pipelines(db, infra)) {
        std::cerr << "Problem detected while loading pipelines";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.sics) ) {
        std::cerr << "Problem detected while loading initial condition for stations";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.pics) ) {
        std::cerr << "Problem detected while loading initial condition for pipelines";
        std::cerr << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

 } // namespace shimmer