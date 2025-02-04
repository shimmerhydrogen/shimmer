-- Insert graph nodes
-- (name, number, height, type) 
-- the type must match the values in station_types
insert into stations(s_number, s_name, t_type) values
    ( 1, 's01_entry', 1),
    ( 2, 's02_exit', 2),
    ( 3, 's03_junct', 3),
    ( 4, 's04_junct', 3),
    ( 5, 's05_junct', 3),
    ( 6, 's06_exit', 2),
    ( 7, 's07_junct', 3),
    ( 8, 's08_junct', 3),
    ( 9, 's09_exit', 2),
    (10, 's10_junct', 3),
    (11, 's11_exit', 2),
    (12, 's12_entry', 1),
    (13, 's13_exit', 2),
    (14, 's14_junct', 3),
    (15, 's15_junct', 3);

select s_name, s_number, t_descr
    from stations inner join station_types
    where stations.t_type = station_types.t_type;

-- Insert graph edges
-- (name, station_from, station_to, type)
-- the type must match the values in pipeline_types
insert into pipelines values
    ('pipe01',  1,  4, 0),
    ('pipe02',  2,  3, 0),
    ('pipe03',  3,  4, 0),
    ('pipe05',  4,  5, 0),
    ('pipe04',  3,  5, 0),
    ('pipe06',  5,  6, 0),
    ('pipe07',  5, 14, 0),
    ('pipe09',  8,  7, 0),
    ('pipe08',  7,  5, 0),
    ('pipe10', 12,  7, 0),
    ('pipe11', 13,  8, 0),
    ('pipe12',  9,  8, 0),
    ('pipe13',  8, 10, 0),
    ('pipe14', 10, 11, 0),
    ('pipe15',  4, 10, 0),
    ('pipe16', 14, 15, 2),
    ('pipe17', 15, 18, 0);

select p_name, s_from, s_to, t_name
    from pipelines inner join pipeline_types
    where pipelines.p_type = pipeline_types.p_type;

insert into pipe_parameters values ('pipe01',  1,  4, 1.2,  80000, 1.20e-5);
insert into pipe_parameters values ('pipe02',  2,  3, 0.6,  16000, 1.20e-5);
insert into pipe_parameters values ('pipe03',  3,  4, 0.8,  40000, 1.20e-5);
insert into pipe_parameters values ('pipe04',  3,  5, 0.7, 160000, 1.20e-5);
insert into pipe_parameters values ('pipe05',  4,  5, 0.8, 200000, 1.20e-5);
insert into pipe_parameters values ('pipe06',  5,  6, 0.6,  24000, 1.20e-5);
insert into pipe_parameters values ('pipe07',  5, 14, 0.2,  60000, 1.20e-5);
insert into pipe_parameters values ('pipe08',  7,  5, 0.9,  80000, 1.20e-5);
insert into pipe_parameters values ('pipe09',  8,  7, 0.7,  64000, 1.20e-5);
insert into pipe_parameters values ('pipe10', 12,  7, 0.6, 240000, 1.20e-5);
insert into pipe_parameters values ('pipe11', 13,  8, 0.2,  28000, 1.20e-5);
insert into pipe_parameters values ('pipe12',  9,  8, 0.9,  80000, 1.20e-5);
insert into pipe_parameters values ('pipe13',  8, 10, 0.7, 160000, 1.20e-5);
insert into pipe_parameters values ('pipe14', 10, 11, 0.3,  40000, 1.20e-5);
insert into pipe_parameters values ('pipe15',  4, 10, 0.9, 320000, 1.20e-5);
insert into pipe_parameters values ('pipe17', 15, 18, 0.2,  60000, 1.20e-5);





-- Profile station 1
insert into STATTYPE1 values ( 1,     0, 7000000 );
insert into STATTYPE1 values ( 1,  3600, 6941666.66666667 );
insert into STATTYPE1 values ( 1,  7200, 6883333.33333333 );
insert into STATTYPE1 values ( 1, 10800, 6825000 );
insert into STATTYPE1 values ( 1, 14400, 6766666.66666667 );
insert into STATTYPE1 values ( 1, 18000, 6708333.33333333 );
insert into STATTYPE1 values ( 1, 21600, 6650000 );
insert into STATTYPE1 values ( 1, 25200, 6708333.33333333 );
insert into STATTYPE1 values ( 1, 28800, 6766666.66666667 );
insert into STATTYPE1 values ( 1, 32400, 6825000 );
insert into STATTYPE1 values ( 1, 36000, 6883333.33333333 );
insert into STATTYPE1 values ( 1, 39600, 6941666.66666667 );
insert into STATTYPE1 values ( 1, 43200, 7000000 );
insert into STATTYPE1 values ( 1, 46800, 7233333.33333333 );
insert into STATTYPE1 values ( 1, 50400, 7466666.66666667 );
insert into STATTYPE1 values ( 1, 54000, 7700000 );
insert into STATTYPE1 values ( 1, 57600, 7933333.33333333 );
insert into STATTYPE1 values ( 1, 61200, 8166666.66666666 );
insert into STATTYPE1 values ( 1, 64800, 8400000 );
insert into STATTYPE1 values ( 1, 68400, 8166666.66666667 );
insert into STATTYPE1 values ( 1, 72000, 7933333.33333333 );
insert into STATTYPE1 values ( 1, 75600, 7700000 );
insert into STATTYPE1 values ( 1, 79200, 7466666.66666667 );
insert into STATTYPE1 values ( 1, 82800, 7233333.33333334 );
insert into STATTYPE1 values ( 1, 86400, 7000000 );

-- Profile station 12
insert into STATTYPE12 values ( 12, 0, 7000000 );
insert into STATTYPE12 values ( 12, 3600, 6766666.66666667 );
insert into STATTYPE12 values ( 12, 7200, 6533333.33333333 );
insert into STATTYPE12 values ( 12, 10800, 6300000 );
insert into STATTYPE12 values ( 12, 14400, 6066666.66666667 );
insert into STATTYPE12 values ( 12, 18000, 5833333.33333333 );
insert into STATTYPE12 values ( 12, 21600, 5600000 );
insert into STATTYPE12 values ( 12, 25200, 5833333.33333333 );
insert into STATTYPE12 values ( 12, 28800, 6066666.66666667 );
insert into STATTYPE12 values ( 12, 32400, 6300000 );
insert into STATTYPE12 values ( 12, 36000, 6533333.33333333 );
insert into STATTYPE12 values ( 12, 39600, 6766666.66666667 );
insert into STATTYPE12 values ( 12, 43200, 7000000 );
insert into STATTYPE12 values ( 12, 46800, 7233333.33333333 );
insert into STATTYPE12 values ( 12, 50400, 7466666.66666667 );
insert into STATTYPE12 values ( 12, 54000, 7700000 );
insert into STATTYPE12 values ( 12, 57600, 7933333.33333333 );
insert into STATTYPE12 values ( 12, 61200, 8166666.66666666 );
insert into STATTYPE12 values ( 12, 64800, 8400000 );
insert into STATTYPE12 values ( 12, 68400, 8458333.33333333 );
insert into STATTYPE12 values ( 12, 72000, 8516666.66666667 );
insert into STATTYPE12 values ( 12, 75600, 8575000 );
insert into STATTYPE12 values ( 12, 79200, 8633333.33333333 );
insert into STATTYPE12 values ( 12, 82800, 8691666.66666667 );
insert into STATTYPE12 values ( 12, 86400, 8750000 );
