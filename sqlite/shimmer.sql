PRAGMA foreign_keys = ON;

create table station_types (
    s_type      INTEGER NOT NULL,
    t_descr     TEXT,
    PRIMARY KEY(s_type)
);

insert into station_types values (0, 'REMI_WO_BACKFLOW');
insert into station_types values (1, 'INJ_W_PRESS_CONTROL');
insert into station_types values (2, 'OUTLET');
insert into station_types values (3, 'JUNCTION');
insert into station_types values (4, 'CONSUMPTION_WO_PRESS');

-- The stations. They are the nodes of the graph
create table stations (
    s_name      TEXT,
    s_number    INTEGER,
    s_height    REAL,
    s_type      INTEGER,
    PRIMARY KEY(s_number),
    
    -- The type of the station must be well-defined
    FOREIGN KEY (s_type)
        REFERENCES station_types(s_type)
);

create table station_fd_parameters (
    s_number    INTEGER,
    s_pressure  REAL,
    s_massflow  REAL,
    
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number)
);


-- Pipeline element type. Can be a pipe, a compressor, a regulator, ...
create table pipeline_types (
    p_type      INTEGER,
    t_name      TEXT NOT NULL,
    PRIMARY KEY (p_type)
);

insert into pipeline_types values (0, 'pipeline');
insert into pipeline_types values (1, 'resistor');
insert into pipeline_types values (2, 'compressor');
insert into pipeline_types values (3, 'regulator');
insert into pipeline_types values (4, 'valve');

-- The pipelines. They are the edges of the graph.
create table pipelines (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,    
    p_type      INTEGER,
    PRIMARY KEY (p_name, s_from, s_to),

    -- The source station must exist
    FOREIGN KEY (s_from)
        REFERENCES stations(s_number),
    -- The destination station must exist
    FOREIGN KEY (s_to)
        REFERENCES stations(s_number)
    -- The pipeline type must be valid
    FOREIGN KEY (p_type)
        REFERENCES pipeline_types(p_type)
);

-- Pipeline parameters as length, diameter and so on.
create table pipeline_parameters (
    p_name      TEXT_NOT_NULL,
    s_from      INTEGER,
    s_to        INTEGER,
    length      REAL NOT NULL,
    diameter    REAL NOT NULL,
    epsi        REAL NOT NULL,

    -- The referenced pipeline must exist
    FOREIGN KEY (p_name, s_from, s_to)
        REFERENCES pipelines(p_name, s_from, s_to)
);

-- The gases. Which are the parameters associated to each gas?
create table gases (
    g_name      TEXT,
    --g_descr     TEXT,
    PRIMARY KEY (g_name)
);

insert into gases(g_name) values ('CH4'), ('N2'), ('CO2'), ('C2H6'), ('C3H8'),
    ('i_C4H10'), ('n_C4H10'), ('i_C5H12'), ('n_C5H12'), ('C6H14'), ('C7H16'),
    ('C8H18'), ('C9H20'), ('C10H22'), ('H2'), ('O2'), ('CO'), ('H2O'), ('H2S'),
    ('He'), ('Ar');


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




insert into stations values ('station 0',   0,  0,  0);
insert into stations values ('station 1',   1,  0,  2);
insert into stations values ('station 2',   2,  0,  3);
insert into stations values ('station 3',   3,  0,  3);
insert into stations values ('station 4',   4,  0,  3);
insert into stations values ('station 5',   5,  0,  4);
insert into stations values ('station 6',   6,  0,  3);
insert into stations values ('station 7',   7,  0,  3);
insert into stations values ('station 8',   8,  0,  4);
insert into stations values ('station 9',   9,  0,  3);
insert into stations values ('station 10', 10,  0,  4);
insert into stations values ('station 11', 11,  0,  1);
insert into stations values ('station 12', 12,  0,  2);

insert into station_fd_parameters values (0, 70.000000000, -75);
insert into station_fd_parameters values (1, 70.000000000,  20);
insert into station_fd_parameters values (2, 69.300000000,   0);
insert into station_fd_parameters values (3, 69.300000000,   0);
insert into station_fd_parameters values (4, 68.607000000,   0);
insert into station_fd_parameters values (5, 67.920930000,  20);
insert into station_fd_parameters values (6, 67.241720700,   0);
insert into station_fd_parameters values (7, 67.920930000,   0);
insert into station_fd_parameters values (8, 67.241720700,  50);
insert into station_fd_parameters values (9, 67.241720700,   0);
insert into station_fd_parameters values (10, 66.569303493,  15);
insert into station_fd_parameters values (11, 70.000000000, -40);
insert into station_fd_parameters values (12, 67.241720700,  10);




select stations.s_name, station_types.t_descr
    from stations
    inner join station_types on stations.s_type = station_types.s_type;

insert into pipelines values ('pipe1', 1, 2, 0);
insert into pipelines values ('pipe2', 1, 3, 0);
insert into pipelines values ('pipe3', 2, 4, 0);

insert into pipeline_parameters values ('pipe1', 1, 2, 100, 10, 0);
insert into pipeline_parameters values ('pipe2', 1, 3, 200, 5, 0);
insert into pipeline_parameters values ('pipe3', 2, 4, 50, 5, 0);

insert into gases values ('gas1'), ('gas2'), ('gas3');

insert into injects values (1, 'gas1', 0.5); 
insert into injects values (1, 'gas2', 0.5); 
insert into injects values (2, 'gas1', 0.3); 
insert into injects values (2, 'gas2', 0.3); 
insert into injects values (2, 'gas3', 0.3); 















-- Ask which gases are injected by station 2
--select stations.s_name, injects.g_name, injects.quantity
--    from stations inner join injects on stations.s_number = injects.s_number
--    where stations.s_number = 2;

-- Ask the total quantity of gas injected by each station
--select stations.s_name, sum(injects.quantity)
--    from stations inner join injects on stations.s_number = injects.s_number
--    group by stations.s_name;

