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

int load(const std::string db_filename, infrastructure& infra)
{
    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
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
    
    if (SHIMMER_SUCCESS != database::load(db, infra.s_u2i, infra.settings_compr_stat)) {
        std::cerr << "Problems detected while loading compressor settings";
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

    if (SHIMMER_SUCCESS != check_pipeline_data_consistency(infra)) {
        std::cerr << "Inconsistencies detected in pipeline data";
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

    //std::cout << "P: " << P.transpose() << std::endl;
    //std::cout << "G: " << G.transpose() << std::endl;
    //std::cout << "L: " << L.transpose() << std::endl;

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

int launch_solver(const config& cfg)
{
    shimmer::infrastructure infra;

    int err = shimmer::load(cfg.database, infra);
    if (err != SHIMMER_SUCCESS) {
        std::cerr << "Problem detected while loading DB" << std::endl;
        return 1;
    }

    shimmer::variable guess = initial_guess(infra);

    /* BEGIN GAS MASS FRACTIONS */
    int nstations = num_stations(infra);
    shimmer::matrix_t y_nodes = shimmer::matrix_t::Zero(nstations, NUM_GASES);
    for (size_t i = 0; i < infra.mass_fractions.size(); i++) {
        const auto& mf = infra.mass_fractions[i];
        assert(mf.i_snum < nstations);
        shimmer::vector_t y = shimmer::vector_t::Zero(NUM_GASES);
        std::copy(mf.fractions.begin(), mf.fractions.end(), y.begin());
        y_nodes.row(i) = y;
    }

    shimmer::incidence inc(infra.graph);
    shimmer::matrix_t y_pipes = inc.matrix_in().transpose() * y_nodes;   
    /* END GAS MASS FRACTIONS */

    using time_solver_t = shimmer::time_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    time_solver_t ts1(infra.graph, cfg.temperature);
    ts1.initialization(guess, cfg.dt_std, cfg.tol_std, y_nodes, y_pipes);  
    ts1.advance(cfg.dt, cfg.steps, cfg.tol, y_nodes, y_pipes);
    auto sol_full  = ts1.solution_full();
    std::cout << sol_full << std::endl;

    std::cout << sol_full.rows() << " " << sol_full.cols() << std::endl;

    auto Pbegin = 0;
    auto Plen = num_stations(infra);
    shimmer::save_pressures(cfg.database, infra,
        sol_full.block(0, Pbegin, sol_full.rows(), Plen)
    );

    auto Lbegin = num_stations(infra);
    auto Llen = num_pipes(infra);
    shimmer::save_flowrates(cfg.database, infra,
        sol_full.block(0, Lbegin, sol_full.rows(), Llen)
    );

    return 0;
}


} // namespace shimmer