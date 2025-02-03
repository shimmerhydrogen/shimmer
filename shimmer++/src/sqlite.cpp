#include <iostream>
#include <cstdio>
#include "boundary.h"

#include "sqlite.hpp"





namespace shimmer {

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
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_));
        sqlite3_close(db_);
        return 1;
    }

    return 0;
}

network_database::~network_database()
{
    sqlite3_close(db_);
}

template<typename T>
concept station_setting = requires {
    T::u_snum;
    T::i_snum;
    T::Pmin;
    T::Pmax;
    T::Lmin;
    T::Lmax;
};

template<station_setting S>
auto
lookup_station_setting(const std::vector<S>& settings, int i_snum)
{
    auto comp = [](const S& s, int i) {
        return s.i_snum < i;
    };

    auto sbegin = settings.begin();
    auto send = settings.end();
    auto sitor = std::lower_bound(sbegin, send, i_snum, comp);
    if ( (sitor == send) or ((*sitor).i_snum != i_snum) )
        return send;

    return sitor;
}

/* Temporary helper */
template<station_setting S>
auto
convert_limits(const S& s)
{
    std::vector<pair_input_t> user_constraints = {
        std::make_pair(P_GREATER_EQUAL, s.Pmin),
        std::make_pair(P_LOWER_EQUAL,   s.Pmax),
        std::make_pair(L_GREATER_EQUAL, s.Lmin),
        std::make_pair(L_LOWER_EQUAL,   s.Lmax)
    };

    return user_constraints;
}

/* Temporary helper */
template<station_setting S>
auto
convert_Lprof(const S& s)
{
    Eigen::VectorXd ret( s.Lprofile.size() );
    for (size_t i = 0; i < s.Lprofile.size(); i++)
        ret[i] = s.Lprofile[i].value;
    return ret;
}

/* Temporary helper */
template<station_setting S>
auto
convert_Pprof(const S& s)
{
    Eigen::VectorXd ret( s.Pprofile.size() );
    for (size_t i = 0; i < s.Pprofile.size(); i++)
        ret[i] = s.Pprofile[i].value;
    return ret;
}

int
network_database::populate_type_dependent_station_data(vertex_properties& vp)
{
    switch (vp.type) {
        
        case station_type_x::JUNCTION: {
            vp.node_station = std::make_unique<junction>();
            break;
        }

        case(station_type_x::REMI_WO_BACKFLOW): {
            auto itor = lookup_station_setting(settings_remi_wo, vp.i_snum);
            if ( itor == settings_remi_wo.end() ) {
                std::cout << "Warning: No data for station " << vp.u_snum;
                std::cout << " (ReMi w/o pressure control)" << std::endl;
                return 1;
            }
            const setting_remi_wo &setting = *itor;
            assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

            auto limits = convert_limits(setting);
            auto Pset = convert_Pprof(setting);
            auto remi = make_remi_wo_backflow(Pset, limits, limits);
            vp.node_station = std::make_unique<multiple_states_station>(remi);                
            break;
        }
            
        case(station_type_x::INJ_W_PRESS_CONTROL): {
            auto itor = lookup_station_setting(settings_injection_w, vp.i_snum);
            if ( itor == settings_injection_w.end() ) {
                std::cout << "Warning: No data for station " << vp.u_snum;
                std::cout << " (Injection w/ pressure control)" << std::endl;
                return 1;
            }
            const setting_injection_w &setting = *itor;
            assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

            auto limits = convert_limits(setting);
            auto Pset = convert_Pprof(setting);
            auto Lset = convert_Lprof(setting);
            auto f = setting.f;
            auto inj_station = make_inj_w_pressure(f, Pset, Lset,
                                              limits,
                                              limits);
            vp.node_station = std::make_unique<multiple_states_station>(inj_station);
            break;
        }

        case(station_type_x::CONSUMPTION_WO_PRESS):
        {
            auto itor = lookup_station_setting(settings_conspoint_wo, vp.i_snum);
            if ( itor == settings_conspoint_wo.end() ) {
                std::cout << "Warning: No data for station " << vp.u_snum;
                std::cout << " (consumption w/o pressure control)" << std::endl;
                return 1;
            }
            const setting_conspoint_wo &setting = *itor;
            assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

            auto limits = convert_limits(setting);
            auto Lset = convert_Lprof(setting);
            auto consumption = make_consumption_wo_press(Lset, limits);
            vp.node_station = std::make_unique<one_state_station>(consumption);
            break;
        }
        case(station_type_x::OUTLET):
        {
            auto itor = lookup_station_setting(settings_outlet, vp.i_snum);
            if ( itor == settings_outlet.end() ) {
                std::cout << "Warning: No data for station " << vp.u_snum;
                std::cout << " (Outlet)" << std::endl;
                return 1;
            }
            const setting_outlet &setting = *itor;
            assert((setting.u_snum == vp.u_snum) and (setting.i_snum == vp.i_snum));

            auto Lset = convert_Lprof(setting);
            auto exit_station = make_outlet(Lset);
            vp.node_station = std::make_unique<one_state_station>(exit_station);
            break;
        }
            
        default:
            throw std::invalid_argument("Invalid station type");  
    }
    return 0;
}

/* Callback for aggregate functions COUNT() and MAX() */
static int
af_callback( void *result_vp, int count, char **data, char **columns) {
    int *result_p = (int *)result_vp;
    *result_p = atoi(data[0]);
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
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }
    
    /* Get the number of stations in the database */
    int station_count;
    q = "SELECT COUNT(s_number) FROM stations";
    rc = sqlite3_exec(db_, q.c_str(), af_callback, &station_count, &zErrMsg);
    if (rc) {
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    /* Resize the mapping arrays */
    s_u2i.resize(max_station_number+1);
    s_i2u.resize(station_count);
    
    /* Compute the actual mapping */
    q = "SELECT * FROM stations ORDER BY s_number";
    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    int i_num = 0;
    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int u_num = sqlite3_column_int(stmt, 0);
        s_u2i.at(u_num) = i_num;
        s_i2u.at(i_num) = u_num;
        i_num++;
    }
    
    sqlite3_finalize(stmt);
    return 0;
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
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    s_u2vd.resize( s_u2i.size() );
    s_i2vd.resize( s_i2u.size() );

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        //int num_cols = sqlite3_column_count(stmt);
        vertex_properties vp;
        
        /* Station number, convert from user to internal */
        int s_un = sqlite3_column_int(stmt, 0);
        int s_in = s_u2i.at(s_un).value();
        vp.node_num = s_in;
        vp.u_snum = s_un;
        vp.i_snum = s_in;

        /* Station name */
        vp.name = (char *) sqlite3_column_text(stmt, 1);
        
        vp.type = static_cast<station_type_x>(sqlite3_column_int(stmt,2));

        /* Station location */
        vp.height = sqlite3_column_double(stmt, 3);
        vp.latitude = sqlite3_column_double(stmt, 4);
        vp.longitude = sqlite3_column_double(stmt, 5);
        
        populate_type_dependent_station_data(vp);
        
        vertex_descriptor vtxd = boost::add_vertex(g);
        /* keep track of the mapping from user/internal numbering to boost
         * vertex descriptors. this should be done with boost::labeled_graph,
         * however it does not like unique pointers, so here we are. */
        s_u2vd.at(s_un) = vtxd;
        s_i2vd.at(s_in) = vtxd;
        
        /* Move the newly created vertex_properties to the graph */
        g[vtxd] = std::move(vp);
    }

    sqlite3_finalize(stmt);

    return 0;
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
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int from_un = sqlite3_column_int(stmt, 1);
        int to_un = sqlite3_column_int(stmt, 2);
        auto from_vtx = s_u2vd.at(from_un).value();
        auto to_vtx = s_u2vd.at(to_un).value();
        
        
        int type = sqlite3_column_int(stmt, 3);
        
        boost::add_edge(from_vtx, to_vtx, g);
    }

    sqlite3_finalize(stmt);

    return 0;
}

int
network_database::populate_graph(infrastructure_graph& g)
{
    if (!db_)
        return 1;

    /* Make the mappings from user numbering to internal numbering */
    renumber_stations();

    /* Import the data for all the stations */
    //import_outlet(settings_outlet);
    import_remi_wo(settings_remi_wo);
    import_injection_w(settings_injection_w);


    /* Import the graph */
    import_stations(g);
    import_pipelines(g);

    return 0;
}

} // namespace shimmer


