create table stations (
    s_name      TEXT NOT NULL,
    s_number    INTEGER,
    s_height    REAL,
    PRIMARY KEY(s_number)
);

create table pipelines (
    p_name      TEXT NOT NULL,
    s_from      INTEGER,
    s_to        INTEGER,    
    lenght      REAL NOT NULL,
    diameter    REAL NOT NULL,
    PRIMARY KEY (p_name, s_from, s_to),
    FOREIGN KEY (s_from)
        REFERENCES stations(s_number),
    FOREIGN KEY (s_to)
        REFERENCES stations(s_number)
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

insert into stations values ('station1', 1, 0);
insert into stations values ('station2', 2, 0);
insert into stations values ('station3', 3, 0);
insert into stations values ('station4', 4, 0);
insert into stations values ('station5', 5, 0);

insert into pipelines values ('pipe1', 1, 2, 100, 10);
insert into pipelines values ('pipe2', 1, 3, 100, 10);
insert into pipelines values ('pipe3', 2, 4, 100, 10);

insert into gases values ('gas1'), ('gas2'), ('gas3');

insert into injects values (1, 'gas1', 0.5); 
insert into injects values (1, 'gas2', 0.5); 
insert into injects values (2, 'gas1', 0.3); 
insert into injects values (2, 'gas2', 0.3); 
insert into injects values (2, 'gas3', 0.3); 

select stations.s_name, injects.g_name, injects.quantity
    from stations inner join injects on stations.s_number = injects.s_number
    where stations.s_number = 2;

select stations.s_name, sum(injects.quantity)
    from stations inner join injects on stations.s_number = injects.s_number
    group by stations.s_name;

