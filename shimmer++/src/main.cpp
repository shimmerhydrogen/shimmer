/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#include <iostream>
#include <fstream>
#include <string>

#include <QApplication>
#include <QtGui>
#include <QWidget>

#include "sol/sol.hpp"

#include "infrastructure_graph.h"
#include "xlsx_io.h"

using namespace shimmer;

static bool
init_lua(sol::state& lua)
{
    lua.open_libraries(sol::lib::base);
    lua.create_table("config");
    sol::load_result script = lua.load_file("../share/shimmer_config.lua");
    if (!script.valid()) {
        std::cout << "Can't load config file" << std::endl;
        return false;
    }
    script();
    return true;
}

static void
make_dummy_infrastructure(infrastructure_graph& igraph)
{
    auto add_vertex = [&](vertex_properties&& vp) 
    {
        auto v = boost::add_vertex(igraph);
        igraph[v] = std::move(vp);

        return v;
    };

    std::vector<vertex_descriptor> vds;
    vds.push_back( add_vertex( vertex_properties( "station1", 1, 100., 10.,0) ) );
    vds.push_back( add_vertex( vertex_properties( "station2", 2, 120., 30.,0) ) );
    vds.push_back( add_vertex( vertex_properties( "station3", 3,  90., 20.,0) ) );
    vds.push_back( add_vertex( vertex_properties( "station4", 4, 160., 40.,0) ) );

    edge_properties ep;
    boost::add_edge( vds[0], vds[1], ep, igraph);
    boost::add_edge( vds[1], vds[2], ep, igraph);
    boost::add_edge( vds[1], vds[3], ep, igraph);
}

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    sol::state lua;
    if (not init_lua(lua))
        return 1;

    infrastructure_graph igraph;

    //import_infrastructure_from_xlsx(lua, igraph);
    make_dummy_infrastructure(igraph);
    write_graphviz("test_graph.dot", igraph);

    return 0;

    QWidget w;
    w.show();
    return app.exec();
}

