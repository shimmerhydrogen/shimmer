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

#include <iostream>
#include <fstream>
#include <string>

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
    sol::state lua;
    if (not init_lua(lua))
        return 1;

    infrastructure_graph igraph;

    //import_infrastructure_from_xlsx(lua, igraph);
    make_dummy_infrastructure(igraph);
    write_graphviz("test_graph.dot", igraph);

    return 0;
}

