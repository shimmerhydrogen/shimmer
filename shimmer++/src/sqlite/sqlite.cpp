#include <iostream>
#include <cstdio>
#include "boundary.h"
#include "errors.h"
#include "sqlite.hpp"

namespace shimmer {

std::optional<table_name_pair_t>
limits_and_profile_table_names(sqlite3 *db, station_type stat_type)
{
    sqlite3_stmt *stmt = nullptr;

    enum class col : int {
        t_limits_table = 0,
        t_profile_table = 1
    };

    std::string q =
        "SELECT t_limits_table, t_profile_table "
        "FROM station_types "
        "WHERE t_type = ?";

    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(db) << std::endl;
        return {};
    }

    rc = sqlite3_bind_int(stmt, 1, +stat_type);

    if ( sqlite3_step(stmt) != SQLITE_ROW ) {
        std::cerr << "Shimmer DB: Invalid station type " << +stat_type << std::endl;
        return {};
    }
    
    std::string limits_table_name;
    if ( sqlite3_column_type(stmt, +col::t_limits_table) == SQLITE3_TEXT ) {
        limits_table_name = (const char *) sqlite3_column_text(stmt, +col::t_limits_table);
    }

    std::string profile_table_name;
    if ( sqlite3_column_type(stmt, +col::t_profile_table) == SQLITE3_TEXT ) {
        profile_table_name = (const char *) sqlite3_column_text(stmt, +col::t_profile_table);
    }

    if ( (limits_table_name == "") or (profile_table_name == "") ) {
        return {};
    }

    rc = sqlite3_clear_bindings(stmt);
    rc = sqlite3_reset(stmt);
    rc = sqlite3_finalize(stmt);

    return std::pair(limits_table_name, profile_table_name);
}

std::optional<std::string>
table_name(sqlite3 *db, setting_table st, station_type stat_type)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q;
    if ( st == setting_table::limits ) {
        q = "SELECT t_limits_table "
        "FROM station_types "
        "WHERE t_type = ?";
    } else if ( st == setting_table::profiles ) {
        q = "SELECT t_profile_table "
        "FROM station_types "
        "WHERE t_type = ?";
    } else {
        std::cerr << "Invalid setting table specified" << std::endl;
        return {};
    }

    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return {};
    }

    rc = sqlite3_bind_int(stmt, 1, +stat_type);

    if ( sqlite3_step(stmt) != SQLITE_ROW ) {
        std::cerr << "Shimmer DB: Invalid station type " << +stat_type << std::endl;
        return {};
    }
    
    std::string tname;
    if ( sqlite3_column_type(stmt, /*column*/ 0) == SQLITE3_TEXT ) {
        tname = (const char *) sqlite3_column_text(stmt, /*column*/ 0);
    }

    if (tname == "") {
        return {};
    }

    rc = sqlite3_clear_bindings(stmt);
    rc = sqlite3_reset(stmt);
    rc = sqlite3_finalize(stmt);

    return tname;
}

network_database::network_database()
    : db_(nullptr)
{}

network_database::network_database(const std::string& filename)
{
    open(filename);
}

int
network_database::open(const std::string& filename)
{
    int rc;
    rc = sqlite3_open(filename.c_str(), &db_);
    if(rc) {
        std::cerr << "Can't open database: " << sqlite3_errmsg(db_) << std::endl;
        sqlite3_close(db_);
        return SHIMMER_DATABASE_PROBLEM;
    }

    return SHIMMER_SUCCESS;
}

network_database::~network_database()
{
    sqlite3_close(db_);
}

template<typename T>
concept limits = requires {
    T::Pmin;
    T::Pmax;
    T::Lmin;
    T::Lmax;
};

/* Temporary helper */
auto
convert_limits(const limits auto& l)
{
    std::vector<pair_input_t> user_constraints = {
        std::make_pair(P_GREATER_EQUAL, l.Pmin),
        std::make_pair(P_LOWER_EQUAL,   l.Pmax),
        std::make_pair(L_GREATER_EQUAL, l.Lmin),
        std::make_pair(L_LOWER_EQUAL,   l.Lmax)
    };

    return user_constraints;
}

template<typename T>
concept profile =
    requires { T::Lprofile; } ||
    requires { T::Pprofile; };

/* Temporary helper */
auto
convert_Lprof(const profile auto& p)
{
    Eigen::VectorXd ret( p.Lprofile.size() );
    for (size_t i = 0; i < p.Lprofile.size(); i++)
        ret[i] = p.Lprofile[i].value;
    return ret;
}

/* Temporary helper */
auto
convert_Pprof(const profile auto& p)
{
    Eigen::VectorXd ret( p.Pprofile.size() );
    for (size_t i = 0; i < p.Pprofile.size(); i++)
        ret[i] = p.Pprofile[i].value;
    return ret;
}

int
network_database::populate_type_dependent_station_data(vertex_properties& vp)
{
    switch (vp.type) {
        
        case station_type::JUNCTION: {
            vp.node_station = std::make_unique<junction>();
            break;
        }

        case(station_type::ENTRY_P_REG): {
            auto itor = lookup(nd_.settings_entry_p_reg, vp.i_snum);
            if ( itor == nd_.settings_entry_p_reg.end() ) {
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
            break;
        }
            
        case(station_type::ENTRY_L_REG): {
            auto itor = lookup(nd_.settings_entry_l_reg, vp.i_snum);
            if ( itor == nd_.settings_entry_l_reg.end() ) {
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
            break;
        }

        case(station_type::EXIT_L_REG): {
            auto itor = lookup(nd_.settings_exit_l_reg, vp.i_snum);
            if ( itor == nd_.settings_exit_l_reg.end() ) {
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
            break;
        }
        
        
        case(station_type::PRIVATE_OUTLET): {
            auto itor = lookup(nd_.settings_outlet, vp.i_snum);
            if ( itor == nd_.settings_outlet.end() ) {
                std::cout << "Warning: No data for station " << vp.u_snum;
                std::cout << " (Outlet)" << std::endl;
                return 1;
            }
            const setting_outlet &setting = *itor;
            assert( setting.i_snum == vp.i_snum );

            auto Lset = convert_Lprof(setting);
            auto exit_station = priv::make_station_outlet(Lset);
            vp.node_station = std::make_unique<one_state_station>(exit_station);
            break;
        }
            
        default:
            std::cerr << "WARNING: Unhandled station type " << +vp.type;
            std::cerr << " while loading type-dependent station data" << std::endl;
            return SHIMMER_INVALID_DATA;
    }
    return SHIMMER_SUCCESS;
}

int
network_database::populate_type_dependent_pipe_data(edge_properties& ep, int i_from, int i_to)
{
    //edge_properties ep15 = {    pipe_type::COMPR_STAT,
    //    15, 1, 0.2 ,	1.20E-05,
    //    };

    //std::unordered_map<external_t, std::pair<edge_constraint_t, double>>  user_limits; 

    //user_limits[external_t::P_OUT_MAX] = std::make_pair(edge_constraint_t::LOWER_EQUAL,   80.E+5);
    //user_limits[external_t::P_IN_MIN]  = std::make_pair(edge_constraint_t::GREATER_EQUAL, 50.E+5);
    //user_limits[external_t::BETA_MAX]  = std::make_pair(edge_constraint_t::LOWER_EQUAL,   2.0);
    //user_limits[external_t::BETA_MIN]  = std::make_pair(edge_constraint_t::GREATER_EQUAL, 1.2);
    //user_limits[external_t::FLUX_MAX]  = std::make_pair(edge_constraint_t::LOWER_EQUAL,   80.0);
    //user_limits[external_t::PWD_NOMINAL] = std::make_pair(edge_constraint_t::LOWER_EQUAL, 10.0);
   
    using edge_constraint_t = edge_station::control::constraint_type;
    using external_t = edge_station::external_type; 
    using user_limits_t = std::unordered_map<external_t, std::pair<edge_constraint_t, double>>;

    switch (ep.type) {
        
        case pipe_type::PIPE: {
            auto sitor = lookup(nd_.settings_pipe, i_from, i_to);
            if (sitor == nd_.settings_pipe.end()) {
                auto u_from = nd_.s_i2u.at(i_from);
                auto u_to = nd_.s_i2u.at(i_to);
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

            auto sitor = lookup(nd_.settings_compr_stat, i_from, i_to);
            if (sitor == nd_.settings_compr_stat.end()) {
                auto u_from = nd_.s_i2u.at(i_from);
                auto u_to = nd_.s_i2u.at(i_to);
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
            auto sitor = lookup(nd_.settings_red_stat, i_from, i_to);
            if (sitor == nd_.settings_red_stat.end()) {
                auto u_from = nd_.s_i2u.at(i_from);
                auto u_to = nd_.s_i2u.at(i_to);
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

std::optional<table_name_pair_t>
network_database::limits_and_profile_table_names(station_type stat_type)
{
    return shimmer::limits_and_profile_table_names(db_, stat_type);
}

/* Callback for aggregate functions COUNT() and MAX() */
static int
af_callback( void *result_vp, int count, char **data, char **columns) {
    int *result_p = (int *)result_vp;
    *result_p = (data[0] != nullptr) ? atoi(data[0]) : 0;
    return 0;
}

/* In the database the numbering of the stations (user numbering) can be
 * non-contiguous, so we want to make it contiguous (internal numbering).
 * The mapping is stored in the two vectors `s_u2i` (user-to-internal) and
 * `s_i2u` (internal-to-user). */
int
network_database::renumber_stations()
{
    char *zErrMsg = nullptr;
    
    /* Get the maximum assigned station number */
    int max_station_number;
    std::string q = "SELECT MAX(s_number) FROM stations";
    int rc = sqlite3_exec(db_, q.c_str(), af_callback, &max_station_number, &zErrMsg);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    /* Get the number of stations in the database */
    int station_count;
    q = "SELECT COUNT(s_number) FROM stations";
    rc = sqlite3_exec(db_, q.c_str(), af_callback, &station_count, &zErrMsg);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    if (station_count == 0) {
        return SHIMMER_SUCCESS;
    }

    /* Resize the mapping arrays */
    nd_.s_u2i.resize(max_station_number+1);
    nd_.s_i2u.resize(station_count);
    nnodes_ = station_count;
    
    /* Compute the actual mapping */
    q = "SELECT * FROM stations ORDER BY s_number";
    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    int i_num = 0;
    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int u_num = sqlite3_column_int(stmt, 0);
        nd_.s_u2i.at(u_num) = i_num;
        nd_.s_i2u.at(i_num) = u_num;
        i_num++;
    }
    
    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;
}

int
network_database::import_stations(infrastructure_graph& g)
{              
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM stations";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << zSql << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    nd_.s_u2vd.resize( nd_.s_u2i.size() );
    nd_.s_i2vd.resize( nd_.s_i2u.size() );

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        //int num_cols = sqlite3_column_count(stmt);
        vertex_properties vp;
        
        /* Station number, convert from user to internal */
        int s_un = sqlite3_column_int(stmt, 0);
        int s_in = nd_.s_u2i.at(s_un).value();
        vp.node_num = s_in;
        vp.u_snum = s_un;
        vp.i_snum = s_in;

        /* Station name */
        vp.name = (char *) sqlite3_column_text(stmt, 1);
        
        vp.type = static_cast<station_type>(sqlite3_column_int(stmt,2));

        if ( s_in >= nd_.mass_fractions.size() ) {
            std::cerr << "No mass fractions for station " << s_un << std::endl;
            return SHIMMER_MISSING_DATA;
        }
        const auto& mfs = nd_.mass_fractions.at(s_in).fractions;
        vp.gas_mixture = vector_t::Zero(mfs.size());
        assert(vp.gas_mixture.size() == NUM_GASES);
        std::copy(mfs.begin(), mfs.end(), vp.gas_mixture.begin());

        /* Station location */
        vp.height = sqlite3_column_double(stmt, 3);
        vp.latitude = sqlite3_column_double(stmt, 4);
        vp.longitude = sqlite3_column_double(stmt, 5);
        
        populate_type_dependent_station_data(vp);
        
        vertex_descriptor vtxd = boost::add_vertex(g);
        /* keep track of the mapping from user/internal numbering to boost
         * vertex descriptors. this should be done with boost::labeled_graph,
         * however it does not like unique pointers, so here we are. */
        nd_.s_u2vd.at(s_un) = vtxd;
        nd_.s_i2vd.at(s_in) = vtxd;
        
        /* Move the newly created vertex_properties to the graph */
        g[vtxd] = std::move(vp);
    }

    sqlite3_finalize(stmt);

    return SHIMMER_SUCCESS;
}

int
network_database::import_pipelines(infrastructure_graph& g)
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM pipelines";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << zSql << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    int branch_num = 0;
    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int from_un = sqlite3_column_int(stmt, 1);
        int to_un = sqlite3_column_int(stmt, 2);
        auto from_vtx = nd_.s_u2vd.at(from_un).value();
        auto to_vtx = nd_.s_u2vd.at(to_un).value();

        auto from_in = nd_.s_u2i.at(from_un).value();
        auto to_in = nd_.s_u2i.at(to_un).value();
        
        /* XXX: NAME CURRENTLY NOT USED - CANNOT HAVE
         *      MULTIPLE PIPES BETWEEN THE SAME TWO NODES
         */
        int type = sqlite3_column_int(stmt, 3);
        
        edge_properties ep;
        ep.type = static_cast<pipe_type>(type);
        ep.branch_num = branch_num++;

        populate_type_dependent_pipe_data(ep, from_in, to_in);

        boost::add_edge(from_vtx, to_vtx, ep, g);
    }
    npipes_ = branch_num;

    sqlite3_finalize(stmt);

    return SHIMMER_SUCCESS;
}

int
network_database::populate_graph(infrastructure_graph& g)
{
    if (!db_) {
        std::cerr << "Shimmer DB: database not opened" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Make the mappings from user numbering to internal numbering */
    renumber_stations();

    /* Import the data for all the stations */
    import_outlet(nd_.settings_outlet);
    import_entry_p_reg(nd_.settings_entry_p_reg);
    import_entry_l_reg(nd_.settings_entry_l_reg);
    import_exit_l_reg(nd_.settings_exit_l_reg);

    import_pipe(nd_.settings_pipe);
    import_compr_stat(nd_.settings_compr_stat);

    import_gas_mass_fractions(nd_.mass_fractions);

    /* Import the graph */
    import_stations(g);
    import_pipelines(g);

    /* Import initial conditions */
    import_station_initial_conditions();
    import_pipe_initial_conditions();

    return SHIMMER_SUCCESS;
}

int
network_database::import_entry_l_reg(std::vector<setting_entry_l_reg>& settings)
{
    return database::load(db_, nd_.s_u2i, settings);
}

int
network_database::import_entry_p_reg(std::vector<setting_entry_p_reg>& settings)
{
    return database::load(db_, nd_.s_u2i, settings);
}

int
network_database::import_exit_l_reg(std::vector<setting_exit_l_reg>& settings)
{
    return database::load(db_, nd_.s_u2i, settings);
}

int
network_database::import_pipe(std::vector<setting_pipe>& settings)
{
    return database::load(db_, nd_.s_u2i, settings);
}

int
network_database::import_station_initial_conditions()
{
    return database::load(db_, nd_.s_u2i, nd_.sics);
}

int
network_database::import_pipe_initial_conditions()
{
    return database::load(db_, nd_.s_u2i, nd_.pics);
}

int
network_database::import_gas_mass_fractions(std::vector<gas_mass_fractions>& fracs)
{
    return database::load(db_, nd_.s_u2i, fracs);
}

int
network_database::import_outlet(std::vector<setting_outlet>& settings)
{
    return database::load(db_, nd_.s_u2i, settings);
}

} // namespace shimmer


