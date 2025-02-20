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
