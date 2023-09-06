create table nodes (
    n_name      TEXT,
    n_number    INTEGER,
    n_height    REAL,
    PRIMARY KEY(n_number)
);

create table edge_types (
    e_type      INTEGER,
    type_name   TEXT NOT NULL,
    PRIMARY KEY (e_type)
);

create table edges (
    e_name      TEXT NOT NULL,
    n_from      INTEGER,
    n_to        INTEGER,    
    e_type      INTEGER,
    PRIMARY KEY (e_name, n_from, n_to),
    FOREIGN KEY (n_from)
        REFERENCES nodes(n_number),
    FOREIGN KEY (n_to)
        REFERENCES nodes(n_number)
    FOREIGN KEY (e_type)
        REFERENCES edge_types(e_type)
);

create table pipeline_parameters (
    e_name      TEXT_NOT_NULL,
    n_from      INTEGER,
    n_to        INTEGER,
    length      REAL NOT NULL,
    diameter    REAL NOT NULL,
    epsi        REAL NOT NULL,
    FOREIGN KEY (e_name, n_from, n_to)
        REFERENCES edges(e_name, n_from, n_to)
);

create table gases (
    g_name      TEXT UNIQUE NOT NULL
);

create table injects (
    s_number    INTEGER,
    g_name      TEXT,
    quantity    REAL,
    FOREIGN KEY (s_number)
        REFERENCES stations(s_number),
    FOREIGN KEY (g_name)
        REFERENCES gases(g_name)
);

insert into edge_types values (0, 'pipeline'), (1, 'resistor'),
    (2, 'compressor'), (3, 'regulator'), (4, 'valve');

insert into nodes values ('station1', 1, 0);
insert into nodes values ('station2', 2, 0);
insert into nodes values ('station3', 3, 0);
insert into nodes values ('station4', 4, 0);
insert into nodes values ('station5', 5, 0);

insert into edges values ('pipe1', 1, 2, 0);
insert into edges values ('pipe2', 1, 3, 0);
insert into edges values ('pipe3', 2, 4, 0);

insert into pipeline_parameters values ('pipe1', 1, 2, 100, 10, 0);
insert into pipeline_parameters values ('pipe2', 1, 3, 200, 5, 0);
insert into pipeline_parameters values ('pipe3', 2, 4, 50, 5, 0);

insert into gases values ('gas1'), ('gas2'), ('gas3');

insert into injects values (1, 'gas1', 0.5); 
insert into injects values (1, 'gas2', 0.5); 
insert into injects values (2, 'gas1', 0.3); 
insert into injects values (2, 'gas2', 0.3); 
insert into injects values (2, 'gas3', 0.3); 

#select stations.s_name, injects.g_name, injects.quantity
#    from stations inner join injects on stations.s_number = injects.s_number
#    where stations.s_number = 2;

#select stations.s_name, sum(injects.quantity)
#    from stations inner join injects on stations.s_number = injects.s_number
#    group by stations.s_name;

