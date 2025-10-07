-- Depends on lua-lsqlite3 and lua-posix
-- On Debian: apt install lua-lsqlite3-dev lua-posix-dev

local shmr = require "shimmer";

local ndfname = "test_net.db";

-- Create a new NDF
shmr.create_ndf( ndfname );

-- Open the NDF just created
local ndf = shmr.open_ndf( ndfname );

-- A limit set that we will apply to EXIT stations
local entry_limits = {};
entry_limits.Lmax = -300.0;
entry_limits.Lmin = -10.0;
entry_limits.Pmin = 6000000.0;
entry_limits.Pmax = 8000000.0;

-- A limit set that we will apply to EXIT stations
local exit_limits = {};
exit_limits.Lmin = 10.0;
exit_limits.Lmax = 300.0;
exit_limits.Pmin = 6000000.0;
exit_limits.Pmax = 8000000.0;

-- A profile for ReMi stations
local remi_profile = {
    {0, 7000000.0}
};

-- A profile for Consumption stations
local cons_profile = {
    {0, 250.0}
};

------------------------------------------------------------------
-- Definition of station 0: number 0, type ReMi, h = 100m
shmr.add_station(ndf, 0, shmr.STATION_REMI, 100.0);
-- Apply the limits defined above
shmr.add_remi_limit(ndf, 0, entry_limits);
-- Apply the profile defined above
shmr.add_remi_profile(ndf, 0, remi_profile);
-- Set the station initial conditions: P and L
shmr.set_sic(ndf, 0, 7000000.0, 200);

------------------------------------------------------------------
-- Definition of station 1: number 1, type Junction, h = 200m
shmr.add_station(ndf, 1, shmr.STATION_JUNC, 200.0);
-- Set the station initial conditions: P and L
shmr.set_sic(ndf, 1, 7000000.0, 200);

------------------------------------------------------------------
-- Definition of station 2: number 2, type Consumption, h = 100m
shmr.add_station(ndf, 2, shmr.STATION_CONS, 100.0);
-- Apply the limits defined above
shmr.add_cons_limit(ndf, 2, exit_limits);
-- Apply the profile defined above
shmr.add_cons_profile(ndf, 2, cons_profile);
-- Set the station initial conditions: P and L
shmr.set_sic(ndf, 2, 7000000.0, 200.0);

------------------------------------------------------------------
-- Add a plain pipe from 0 to 1
shmr.add_pipe(ndf, 0, 1);
-- Set params for pipe from 0 to 1: diam = 0.9644m, len = 75000m,
-- roughness = 1.4e-5
shmr.set_pipe_params(ndf, 0, 1, 0.9664, 75000, 1.4e-5);
-- Set initial conditions for pipe from 0 to 1: G
shmr.set_pic(ndf, 0, 1, 200.0);

-- Add a plain pipe from 1 to 2
shmr.add_pipe(ndf, 1, 2);
-- Set params for pipe from 1 to 2: diam = 0.9644m, len = 75000m,
-- roughness = 1.4e-5
shmr.set_pipe_params(ndf, 1, 2, 0.9664, 75000, 1.4e-5);
-- Set initial conditions for pipe from 1 to 2: G
shmr.set_pic(ndf, 1, 2, 200.0);

ndf:close();
