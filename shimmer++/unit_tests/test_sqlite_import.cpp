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
#include "infrastructure_graph.h"
#include "sqlite/sqlite.hpp"

using namespace shimmer;

int main(int argc, char ** argv) {

    if( argc != 2 ) {
        // To make the db: sqlite3 shimmer.db < ../../../sqlite/shimmer.sql
        fprintf(stderr, "Usage: %s DATABASE\n", argv[0]);
        return 1;
    }
  
    network_database db(argv[1]);
    infrastructure_graph g;
    db.populate_graph(g);

    return 0;
}
