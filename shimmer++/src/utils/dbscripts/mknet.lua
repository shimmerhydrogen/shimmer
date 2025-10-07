local shmr = require "shimmer";

local ndfname = "test_net.db";

shmr.create_ndf( ndfname );

local ndf = shmr.open_ndf( ndfname );

local entry_limits = {};
entry_limits.Lmax = -300.0;
entry_limits.Lmin = -10.0;
entry_limits.Pmin = 6000000.0;
entry_limits.Pmax = 8000000.0;

local exit_limits = {};
exit_limits.Lmin = 10.0;
exit_limits.Lmax = 300.0;
exit_limits.Pmin = 6000000.0;
exit_limits.Pmax = 8000000.0;

local remi_profile = {
    {0, 7000000.0}
};

local cons_profile = {
    {0, 250.0}
};

-- Station 0
shmr.add_station(ndf, 0, shmr.STATION_REMI, 100.0);
shmr.add_remi_limit(ndf, 0, entry_limits);
shmr.add_remi_profile(ndf, 0, remi_profile);
shmr.set_sic(ndf, 0, 7000000.0, 200);

-- Station 1
shmr.add_station(ndf, 1, shmr.STATION_JUNC, 200.0);
shmr.set_sic(ndf, 1, 7000000.0, 200);

-- Station 2
shmr.add_station(ndf, 2, shmr.STATION_CONS, 100.0);
shmr.add_cons_limit(ndf, 2, exit_limits);
shmr.add_cons_profile(ndf, 2, cons_profile);
shmr.set_sic(ndf, 2, 7000000.0, 200.0);

-- Pipe from 0 to 1
shmr.add_pipe(ndf, 0, 1);
shmr.set_pipe_params(ndf, 0, 1, 0.9664, 75000, 1.4e-5);
shmr.set_pic(ndf, 0, 1, 200.0);

-- Pipe from 1 to 2
shmr.add_pipe(ndf, 1, 2);
shmr.set_pipe_params(ndf, 1, 2, 0.9664, 75000, 1.4e-5);
shmr.set_pic(ndf, 1, 2, 200.0);

ndf:close();
