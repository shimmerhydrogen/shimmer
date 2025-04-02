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
    (10, 'Inlet station', NULL, NULL),
    (11, 'Outlet station', NULL, NULL);

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

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);

-- Profiles
--  s_number:   number of the station
--  prf_time:   relative time of the sample
--  prf_Pset:   pressure setpoint at the specified time

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
-- Pipeline element type. Can be a pipe, a compressor, a regulator, ...
create table pipeline_types (
    p_type      INTEGER,
    t_name      TEXT NOT NULL,
    PRIMARY KEY (p_type)
);

insert into pipeline_types values
    (0, 'Plain pipe')
    (1, 'Compressor');
    (2, 'Reduction station');
    (3, 'Valve');

-- The pipelines. They are the edges of the graph.
create table pipelines (
    p_name      TEXT NOT NULL,
    s_from      INTEGER NOT NULL,
    s_to        INTEGER NOT NULL,    
    p_type      INTEGER,
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

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
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
create table station_initial_conditions (
    s_number    INTEGER UNIQUE NOT NULL,
    init_P      REAL DEFAULT 0.0 NOT NULL,
    init_L      REAL DEFAULT 0.0 NOT NULL,

    FOREIGN KEY (s_number)
        REFERENCES stations(s_number),
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


-- The gases. Which are the parameters associated to each gas?
create table gases (
    g_name      TEXT,
    --g_descr     TEXT,
    PRIMARY KEY (g_name)
);

insert into gases(g_name) values ('CH4'), ('N2'), ('CO2'), ('C2H6'),
    ('C3H8'), ('i_C4H10'), ('n_C4H10'), ('i_C5H12'), ('n_C5H12'),
    ('C6H14'), ('C7H16'), ('C8H18'), ('C9H20'), ('C10H22'), ('H2'),
    ('O2'), ('CO'), ('H2O'), ('H2S'), ('He'), ('Ar');


-- Who injects what
create table injects (
    s_number    INTEGER,
    g_name      TEXT,
    quantity    REAL,

    -- station number must be valid
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number),
    -- gas name must be valid
    FOREIGN KEY (g_name)
        REFERENCES gases(g_name)
);


--insert into station_parameters values (0, 70.000000000, -75);
--insert into station_parameters values (1, 70.000000000,  20);
--insert into station_parameters values (2, 69.300000000,   0);
--insert into station_parameters values (3, 69.300000000,   0);
--insert into station_parameters values (4, 68.607000000,   0);
--insert into station_parameters values (5, 67.920930000,  20);
--insert into station_parameters values (6, 67.241720700,   0);
--insert into station_parameters values (7, 67.920930000,   0);
--insert into station_parameters values (8, 67.241720700,  50);
--insert into station_parameters values (9, 67.241720700,   0);
--insert into station_parameters values (10, 66.569303493,  15);
--insert into station_parameters values (11, 70.000000000, -40);
--insert into station_parameters values (12, 67.241720700,  10);




select stations.s_name, station_types.t_descr
    from stations
    inner join station_types on stations.t_type = station_types.t_type;

--insert into pipelines values ('pipe1', 1, 2, 0);
--insert into pipelines values ('pipe2', 1, 3, 0);
--insert into pipelines values ('pipe3', 2, 4, 0);

--insert into pipeline_parameters values ('pipe1', 1, 2, 100, 10, 0);
--insert into pipeline_parameters values ('pipe2', 1, 3, 200, 5, 0);
--insert into pipeline_parameters values ('pipe3', 2, 4, 50, 5, 0);

insert into gases values ('gas1'), ('gas2'), ('gas3');

--insert into injects values (1, 'gas1', 0.5); 
--insert into injects values (1, 'gas2', 0.5); 
--insert into injects values (2, 'gas1', 0.3); 
--insert into injects values (2, 'gas2', 0.3); 
--insert into injects values (2, 'gas3', 0.3); 
















-- Ask which gases are injected by station 2
--select stations.s_name, injects.g_name, injects.quantity
--    from stations inner join injects on stations.s_number = injects.s_number
--    where stations.s_number = 2;

-- Ask the total quantity of gas injected by each station
--select stations.s_name, sum(injects.quantity)
--    from stations inner join injects on stations.s_number = injects.s_number
--    group by stations.s_name;

