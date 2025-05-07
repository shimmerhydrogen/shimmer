/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 *
 * The authors (C) 2023, 2024, 2025
 */

#pragma once
#include <concepts>

namespace shimmer {

enum class station_type : int {
    ENTRY_P_REG = 1,        /* ReMi w/o backflow */
    ENTRY_L_REG = 2,        /* Injection w/ pressure control */
    EXIT_L_REG  = 3,        /* Consumption point w/o pressure control */
    JUNCTION    = 4,        /* Junction */
    PRIVATE_INLET  = 10,    /* Inlet, internal use only */
    PRIVATE_OUTLET = 11     /* Outlet, internal use only */
};

enum class pipe_type : int {
    PIPE        = 0,    /* Plain pipe */
    COMPR_STAT  = 1,    /* Compressor station */
    RED_STAT    = 2,    /* Reduction station */
    VALVE       = 3     /* Valve */
};

enum class compressor_mode : int {
    ON_POWER    = 0,    /* Compressor on, Control mode power driver */
    ON_OPRESS   = 1,    /* Compressor on, Control mode outlet pressure */
    ON_IPRESS   = 2,    /* Compressor on, Control mode inlet pressure */
    ON_RATIO    = 3,    /* Compressor on, Control mode compression ratio */
    ON_MASSFLOW = 4,    /* Compressor on, Control mode mass flow */
    OFF_BYPASS  = 10,   /* Compressor off, bypass */
    OFF_CLOSED  = 11    /* Compressor off, closed */
};

enum class reduction_mode : int {
    ON_OPRESS   = 0,    /* Reduction station on, Control mode outlet pressure */
    ON_IPRESS   = 1,    /* Reduction station on, Control mode inlet pressure */
    ON_RATIO    = 2,    /* Reduction station on, Control mode compression ratio */
    ON_MASSFLOW = 3,    /* Reduction station on, Control mode mass flow */
    OFF_BYPASS  = 10,   /* Reduction station off, bypass */
    OFF_CLOSED  = 11    /* Reduction station off, closed */
};

} // namespace shimmer

template<typename T>
    requires std::is_enum_v<T>
constexpr auto operator+(T e) {
    return std::underlying_type_t<T>(e);
}
