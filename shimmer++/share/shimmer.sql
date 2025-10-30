-- SHIMMER database schema.

PRAGMA foreign_keys = ON;

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Descriptions for the types of stations we handle.
--  t_type:     numeric value for the type
--  t_descr:    description for this type of station

create table station_types (
    t_type          INTEGER,
    t_descr         TEXT NOT NULL,
    t_limits_table  TEXT,
    t_profile_table TEXT,

    PRIMARY KEY(t_type)
);

-- Populate with the known station types
insert into station_types values
    (1, 'ReMi station w/o backflow', 'limits_remi_wo', 'profiles_remi_wo'),
    (2, 'Injection station w/ pressure control', 'limits_injection_w', 'profiles_injection_w'),
    (3, 'Consumption point w/o pressure control', 'limits_conspoint_wo', 'profiles_conspoint_wo'),
    (4, 'Junction', NULL, NULL),
    (10, 'Inlet station - private', NULL, NULL),
    (11, 'Outlet station - private', 'limits_outlet_priv', 'profiles_outlet_priv');

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- The stations, or the nodes of the graph.
--  s_name:     name of the station
--  s_number:   number identifying the station
--  s_height:   altimetric height of the station
--  t_type:     type of the station

create table stations (
    s_number    INTEGER,
    s_name      TEXT NOT NULL,
    t_type      INTEGER,

    s_height    REAL DEFAULT 0.0 NOT NULL,
    s_latitude  REAL DEFAULT 0.0 NOT NULL,
    s_longitude REAL DEFAULT 0.0 NOT NULL,

    PRIMARY KEY(s_number),
    
    -- The type of the station must be well-defined
    FOREIGN KEY (t_type)
        REFERENCES station_types(t_type),
    
    CHECK(s_number >= 0)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Limits and profiles.
-- This data for now is more or less the same for all types of
-- station, however we want to keep it in different tables in order to
-- limit the propagation of modifications if in future the data model
-- changes. This is reflected on the C++ side: each type of station
-- has its own I/O code which is roughly the same.

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Limits and profiles for ENTRY stations of type "ReMi without
-- backflow". Parameters described in "Nodal BC.pdf", ENTRY station
-- section 1.1.
--
-- Limits
--  s_number:   number of the station
--  lim_Lmin:   minimum allowed mass flow rate
--  lim_Lmax:   maximum allowed mass flow rate
--  lim_Pmin:   minimum allowed pressure
--  lim_Pmax:   maximum allowed pressure

create table limits_remi_wo (
    s_number    INTEGER UNIQUE,
    lim_Lmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Lmax    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmax    REAL DEFAULT 0.0 NOT NULL,

    CHECK(lim_Lmin <= 0),
    CHECK(lim_Lmax <= lim_Lmin),
    CHECK(lim_Pmin <= lim_Pmax)

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-- Profiles
--  s_number:   number of the station
--  prf_time:   relative time of the sample
--  prf_Pset:   pressure setpoint at the specified time

create table profiles_remi_wo (
    s_number    INTEGER,
    prf_time    REAL DEFAULT 0.0 NOT NULL,
    prf_Pset    REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Limits and profiles for ENTRY stations of type "Injection with
-- pressure control". Parameters described in "Nodal BC.pdf", ENTRY
-- station section 2.2.
--
-- Limits
--  s_number:   number of the station
--  lim_Lmin:   minimum allowed mass flow rate
--  lim_Lmax:   maximum allowed mass flow rate
--  lim_Pmin:   minimum allowed pressure
--  lim_Pmax:   maximum allowed pressure
--  parm_f:     scale factor named "f" in the slides

create table limits_injection_w (
    s_number    INTEGER UNIQUE,
    lim_Lmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Lmax    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmax    REAL DEFAULT 0.0 NOT NULL,
    parm_f      REAL DEFAULT 1.0 NOT NULL,

    CHECK(lim_Lmin <= 0),
    CHECK(lim_Lmax <= lim_Lmin),
    CHECK(lim_Pmin <= lim_Pmax)

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-- Profiles
--  s_number:   number of the station
--  prf_time:   relative time of the sample
--  prf_Lset:   mass flow rate setpoint at the specified time

create table profiles_injection_w (
    s_number    INTEGER,
    prf_time    REAL DEFAULT 0.0 NOT NULL,
    prf_Pset    REAL DEFAULT 0.0 NOT NULL,
    prf_Lset    REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Limits for  (Slides 2.1)
-- Limits and profiles for EXIT stations of type "consumption points
-- w/o pressure control". Parameters described in "Nodal BC.pdf", EXIT
-- station section 1.2.
--
-- Limits
--  s_number:   number of the station
--  lim_Lmin:   minimum allowed mass flow rate
--  lim_Lmax:   maximum allowed mass flow rate
--  lim_Pmin:   minimum allowed pressure
--  lim_Pmax:   maximum allowed pressure

create table limits_conspoint_wo (
    s_number    INTEGER UNIQUE,
    lim_Lmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Lmax    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmin    REAL DEFAULT 0.0 NOT NULL,
    lim_Pmax    REAL DEFAULT 0.0 NOT NULL,

    CHECK(lim_Lmin >= 0),
    CHECK(lim_Lmin <= lim_Lmax),
    CHECK(lim_Pmin <= lim_Pmax),

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-- Profiles
--  s_number:   number of the station
--  prf_time:   relative time of the sample
--  prf_Lset:   pressure setpoint at the specified time

create table profiles_conspoint_wo (
    s_number    INTEGER,
    prf_time    REAL DEFAULT 0.0 NOT NULL,
    prf_Lset    REAL DEFAULT 0.0 NOT NULL,

    CHECK(prf_Lset >= 0),

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Outlet station. This is private, it will disappear in a future release

create table profiles_outlet_priv (
    s_number    INTEGER,
    prf_time    REAL DEFAULT 0.0 NOT NULL,
    prf_Lset    REAL DEFAULT 0.0 NOT NULL,

    CHECK(prf_Lset >= 0),

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Pipeline element type. Can be a pipe, a compressor, a regulator, ...
create table pipeline_types (
    p_type      INTEGER,
    t_name      TEXT NOT NULL,
    PRIMARY KEY (p_type)
);

insert into pipeline_types values
    (0, 'Plain pipe'),
    (1, 'Compressor'),
    (2, 'Reduction station'),
    (3, 'Valve');

-- The pipelines. They are the edges of the graph.
create table pipelines (
    p_name      TEXT NOT NULL,
    s_from      INTEGER NOT NULL,
    s_to        INTEGER NOT NULL,    
    p_type      INTEGER NOT NULL,
    PRIMARY KEY (p_name, s_from, s_to),

    -- The source station must exist
    FOREIGN KEY (s_from)
        REFERENCES stations(s_number),
    -- The destination station must exist
    FOREIGN KEY (s_to)
        REFERENCES stations(s_number),
    -- The pipeline type must be valid
    FOREIGN KEY (p_type)
        REFERENCES pipeline_types(p_type)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Pipe parameters as length, diameter and so on.
create table pipe_parameters (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,
    diameter    REAL DEFAULT 0.0,
    length      REAL DEFAULT 0.0,
    roughness   REAL DEFAULT 0.0,
    ref_nsegs   INTEGER DEFAULT 0,

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to),

    CHECK(ref_nsegs >= 0)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Compressor parameters.
create table compressor_profile (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,

    prf_time    REAL DEFAULT 0.0,
    controlmode INTEGER DEFAULT 10, -- default OFF BYPASS
    power       REAL DEFAULT 0.0,
    outpress    REAL DEFAULT 0.0,
    inpress     REAL DEFAULT 0.0,
    ratio       REAL DEFAULT 0.0,
    massflow    REAL DEFAULT 0.0,

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);

create table compressor_limits (
    p_name          TEXT NOT NULL,
    s_from          INTEGER,
    s_to            INTEGER,

    max_power       REAL DEFAULT 0.0,
    max_outpress    REAL DEFAULT 0.0,
    min_inpress     REAL DEFAULT 0.0,
    max_ratio       REAL DEFAULT 0.0,
    min_ratio       REAL DEFAULT 0.0,
    max_massflow    REAL DEFAULT 0.0,

    PRIMARY KEY (p_name, s_from, s_to),

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);

-----------------------------------------------------------------------
-----------------------------------------------------------------------
-----------------------------------------------------------------------
-----------------------------------------------------------------------
-- Valve parameters.
create table valve_profile (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,

    prf_time    REAL DEFAULT 0.0,
    controlmode INTEGER DEFAULT 10, -- default OFF BYPASS

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);


-----------------------------------------------------------------------
-----------------------------------------------------------------------
create table station_initial_conditions (
    s_number    INTEGER UNIQUE NOT NULL,
    init_P      REAL DEFAULT 0.0 NOT NULL,
    init_L      REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

create table pipe_initial_conditions (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,
    init_G      REAL DEFAULT 0.0 NOT NULL,
    
    PRIMARY KEY (p_name, s_from, s_to),

    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);


create table gas_molar_fractions (
    s_number        INTEGER UNIQUE NOT NULL,
    frac_CH4        REAL DEFAULT 0.0 NOT NULL,
    frac_N2         REAL DEFAULT 0.0 NOT NULL,
    frac_CO2        REAL DEFAULT 0.0 NOT NULL,
    frac_C2H6       REAL DEFAULT 0.0 NOT NULL,
    frac_C3H8       REAL DEFAULT 0.0 NOT NULL,
    frac_iC4H10     REAL DEFAULT 0.0 NOT NULL,
    frac_nC4H10     REAL DEFAULT 0.0 NOT NULL,
    frac_iC5H12     REAL DEFAULT 0.0 NOT NULL,
    frac_nC5H12     REAL DEFAULT 0.0 NOT NULL,
    frac_C6H14      REAL DEFAULT 0.0 NOT NULL,
    frac_C7H16      REAL DEFAULT 0.0 NOT NULL,
    frac_C8H18      REAL DEFAULT 0.0 NOT NULL,
    frac_C9H20      REAL DEFAULT 0.0 NOT NULL,
    frac_C10H22     REAL DEFAULT 0.0 NOT NULL,
    frac_H2         REAL DEFAULT 0.0 NOT NULL,
    frac_O2         REAL DEFAULT 0.0 NOT NULL,
    frac_CO         REAL DEFAULT 0.0 NOT NULL,
    frac_H2O        REAL DEFAULT 0.0 NOT NULL,
    frac_H2S        REAL DEFAULT 0.0 NOT NULL,
    frac_He         REAL DEFAULT 0.0 NOT NULL,
    frac_Ar         REAL DEFAULT 0.0 NOT NULL,
    
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);



create table solution_station_pressures (
    s_number        INTEGER NOT NULL,
    timestep        INTEGER NOT NULL,
    pressure        REAL DEFAULT 0.0 NOT NULL,
    
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

create table solution_pipe_flowrates (
    p_name      TEXT NOT NULL,
    s_from      INTEGER NOT NULL,
    s_to        INTEGER NOT NULL,
    timestep    INTEGER NOT NULL,
    flowrate    REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);

create table solution_station_flowrates (
    s_number        INTEGER NOT NULL,
    timestep        INTEGER NOT NULL,
    flowrate        REAL DEFAULT 0.0 NOT NULL,
    
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

create table solution_pipe_velocities (
    p_name      TEXT NOT NULL,
    s_from      INTEGER NOT NULL,
    s_to        INTEGER NOT NULL,
    timestep    INTEGER NOT NULL,
    velocity    REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);


create table solution_station_molarfrac (
    s_number        INTEGER NOT NULL,
    timestep        INTEGER NOT NULL,
    g_name          TEXT,
    molarfrac       REAL DEFAULT 0.0 NOT NULL,
    
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);


create table gases (
    g_num       INTEGER,
    g_formula   TEXT NOT NULL,
    g_name      TEXT NOT NULL,

    PRIMARY KEY (g_num)

);

insert into gases values
    ( 0, 'CH4',     'Methane'),
    ( 1, 'N2',      'Nitrogen'),
    ( 2, 'CO2',     'Carbon dioxide'),
    ( 3, 'C2H6',    'Ethane'),
    ( 4, 'C3H8',    'Propane'),
    ( 5, 'i_C4H10', 'i-butane'),
    ( 6, 'n_C4H10', 'n-butane'),
    ( 7, 'i_C5H12', 'i-pentane'),
    ( 8, 'n_C5H12', 'n-pentane'),
    ( 9, 'C6H14',   'Hexane'),
    (10, 'C7H16',   'Heptane'),
    (11, 'C8H18',   'Octane'),
    (12, 'C9H20',   'Nonane'),
    (13, 'C10H22',  'Decane'),
    (14, 'H2',      'Hydrogen'),
    (15, 'O2',       'Oxygen'),
    (16, 'CO',      'Carbon oxide'),
    (17, 'H2O',     'Water'),
    (18, 'H2S',      'Hydrogen sulfide'),
    (19, 'He',      'Helium'),
    (20, 'Ar',      'Argon');

