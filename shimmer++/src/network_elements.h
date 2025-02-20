/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 *
 * The authors (C) 2023, 2024, 2025
 */

#pragma once

#if 0
enum station_type {
    ENTRY_P_REG = 1,    /* ReMi w/o backflow */
    ENTRY_L_REG = 2,    /* Injection w/ pressure control */
    EXIT_L_REG  = 3     /* Consumption point w/o pressure control */
};
#endif

enum pipe_type {

};

/* this must become pipe_type */
enum class edge_type {
    pipe,
    compressor,
    regulator,
    valve,
};