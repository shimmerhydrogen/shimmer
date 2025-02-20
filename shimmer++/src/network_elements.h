/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 *
 * The authors (C) 2023, 2024, 2025
 */

#pragma once

#if 0
enum station_type {
    entry_p_reg = 1,    /* ReMi w/o backflow */
    entry_l_reg = 2,    /* Injection w/ pressure control */
    exit_l_reg = 3      /* Consumption point w/o pressure control */
};
#endif

enum pipe_type {

};

enum class edge_type {
    pipe,
    compressor,
    regulator,
    valve,
};