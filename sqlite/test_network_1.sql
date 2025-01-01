-- Insert graph nodes
-- (name, number, height, type) 
-- the type must match the values in station_types
insert into stations(s_number, s_name, t_type) values
    ( 0, 'station 0', 0),
    ( 1, 'station 1', 2),
    ( 2, 'station 2', 3),
    ( 3, 'station 3', 3),
    ( 4, 'station 4', 3),
    ( 5, 'station 5', 4),
    ( 6, 'station 6', 3),
    ( 7, 'station 7', 3),
    ( 8, 'station 8', 4),
    ( 9, 'station 9', 3),
    (10, 'station 10', 4),
    (11, 'station 11', 1),
    (12, 'station 12', 2);

--select s_name, s_number, t_descr
--    from stations inner join station_types
--    where stations.t_type = station_types.t_type;

-- Insert graph edges
-- (name, station_from, station_to, type)
-- the type must match the values in pipeline_types
insert into pipelines values ('pipe-0-3', 0, 3, 0);
insert into pipelines values ('pipe-1-2', 1, 2, 0);
insert into pipelines values ('pipe-2-3', 2, 3, 0);
insert into pipelines values ('pipe-2-4', 2, 4, 0);
insert into pipelines values ('pipe-3-4', 3, 4, 0);
insert into pipelines values ('pipe-4-5', 4, 5, 0);
insert into pipelines values ('pipe-4-7', 4, 7, 0);
insert into pipelines values ('pipe-6-4', 6, 4, 0);
insert into pipelines values ('pipe-7-6', 7, 6, 0);
insert into pipelines values ('pipe-11-6', 11, 6, 0);
insert into pipelines values ('pipe-12-7', 12, 7, 0);
insert into pipelines values ('pipe-8-7', 8, 7, 0);
insert into pipelines values ('pipe-7-9', 7, 9, 0);
insert into pipelines values ('pipe-9-10', 9, 10, 0);
insert into pipelines values ('pipe-3-9', 3, 9, 0);

--select p_name, s_from, s_to from pipelines;

-----------------------------------------------------
-- Station 0 settings

insert into limits_remi_wo values
    (0, -300, -10, 60e5, 80e5);

insert into profiles_remi_wo values
    (0,     0, 7.000000000000000e+06),
    (0,  3600, 6.941666666666667e+06),
    (0,  7200, 6.883333333333334e+06),
    (0, 10800, 6.825000000000000e+06),
    (0, 14400, 6.766666666666667e+06),
    (0, 18000, 6.708333333333334e+06),
    (0, 21600, 6.650000000000000e+06),
    (0, 25200, 6.708333333333333e+06),
    (0, 28800, 6.766666666666666e+06),
    (0, 32400, 6.824999999999998e+06),
    (0, 36000, 6.883333333333333e+06),
    (0, 39600, 6.941666666666666e+06),
    (0, 43200, 7.000000000000000e+06),
    (0, 46800, 7.233333333333333e+06),
    (0, 50400, 7.466666666666666e+06),
    (0, 54000, 7.699999999999997e+06),
    (0, 57600, 7.933333333333330e+06),
    (0, 61200, 8.166666666666663e+06),
    (0, 64800, 8.400000000000000e+06),
    (0, 68400, 8.166666666666666e+06),
    (0, 72000, 7.933333333333333e+06),
    (0, 79200, 7.700000000000000e+06),
    (0, 82800, 7.466666666666669e+06),
    (0, 86400, 7.233333333333336e+06),
    (0, 90000, 7.000000000000000e+06);

-----------------------------------------------------
-- Station 11 settings

insert into limits_injection_w values
    (11, -300, -10, 60e5, 80e5, 1);

insert into profiles_injection_w values
    (11,     0, 7.000000000000000e+06),
    (11,  3600, 6.941666666666667e+06),
    (11,  7200, 6.883333333333334e+06),
    (11, 10800, 6.825000000000000e+06),
    (11, 14400, 6.766666666666667e+06),
    (11, 18000, 6.708333333333334e+06),
    (11, 21600, 6.650000000000000e+06),
    (11, 25200, 6.708333333333333e+06),
    (11, 28800, 6.766666666666666e+06),
    (11, 32400, 6.824999999999998e+06),
    (11, 36000, 6.883333333333333e+06),
    (11, 39600, 6.941666666666666e+06),
    (11, 43200, 7.000000000000000e+06),
    (11, 46800, 7.233333333333333e+06),
    (11, 50400, 7.466666666666666e+06),
    (11, 54000, 7.699999999999997e+06),
    (11, 57600, 7.933333333333330e+06),
    (11, 61200, 8.166666666666663e+06),
    (11, 64800, 8.400000000000000e+06),
    (11, 68400, 8.166666666666666e+06),
    (11, 72000, 7.933333333333333e+06),
    (11, 79200, 7.700000000000000e+06),
    (11, 82800, 7.466666666666669e+06),
    (11, 86400, 7.233333333333336e+06),
    (11, 90000, 7.000000000000000e+06);













