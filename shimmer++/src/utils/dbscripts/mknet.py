# This is the SHIMMER gas network simulator.
# Copyright (C) 2023-2024-2025 Politecnico di Torino
# 
# Dipartimento di Matematica "G. L. Lagrange" - DISMA
# Dipartimento di Energia "G. Ferraris" - DENERG
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import shimmer as shmr

ndfname = "test_net.db"

# Create a new NDF
shmr.create_ndf(ndfname)

# Open the NDF just created
ndf = shmr.open_ndf(ndfname)

# A limit set that we will apply to ENTRY stations
entry_limits = {}
entry_limits['Lmax'] = -300.0
entry_limits['Lmin'] = -10.0
entry_limits['Pmin'] = 6000000.0
entry_limits['Pmax'] = 8000000.0

# A limit set that we will apply to EXIT stations
exit_limits = {}
exit_limits['Lmin'] = 10.0
exit_limits['Lmax'] = 300.0
exit_limits['Pmin'] = 6000000.0
exit_limits['Pmax'] = 8000000.0

# A profile for ReMi stations
remi_profile = [
    (0, 7000000.0)
]

# A profile for Consumption stations
cons_profile = [
    (0, 250.0)
]

# ------------------------------------------------------------------
# Definition of station 0: number 0, type ReMi, h = 100m
shmr.add_station(ndf, 0, shmr.STATION_REMI, 100.0)
# Apply the limits defined above
shmr.add_remi_limit(ndf, 0, entry_limits)
# Apply the profile defined above
shmr.add_remi_profile(ndf, 0, remi_profile)
# Set the station initial conditions: P and L
shmr.set_sic(ndf, 0, 7000000.0, 200)

# ------------------------------------------------------------------
# Definition of station 1: number 1, type Junction, h = 200m
shmr.add_station(ndf, 1, shmr.STATION_JUNC, 200.0)
# Set the station initial conditions: P and L
shmr.set_sic(ndf, 1, 7000000.0, 200)

# ------------------------------------------------------------------
# Definition of station 2: number 2, type Consumption, h = 100m
shmr.add_station(ndf, 2, shmr.STATION_CONS, 100.0)
# Apply the limits defined above
shmr.add_cons_limit(ndf, 2, exit_limits)
# Apply the profile defined above
shmr.add_cons_profile(ndf, 2, cons_profile)
# Set the station initial conditions: P and L
shmr.set_sic(ndf, 2, 7000000.0, 200.0)

# ------------------------------------------------------------------
# Add a plain pipe from 0 to 1
shmr.add_pipe(ndf, 0, 1)
# Set params for pipe from 0 to 1: diam = 0.9644m, len = 75000m,
# roughness = 1.4e-5
shmr.set_pipe_params(ndf, 0, 1, 0.9664, 75000, 1.4e-5)
# Set initial conditions for pipe from 0 to 1: G
shmr.set_pic(ndf, 0, 1, 200.0)

# Add a plain pipe from 1 to 2
shmr.add_pipe(ndf, 1, 2)
# Set params for pipe from 1 to 2: diam = 0.9644m, len = 75000m,
# roughness = 1.4e-5
shmr.set_pipe_params(ndf, 1, 2, 0.9664, 75000, 1.4e-5)
# Set initial conditions for pipe from 1 to 2: G
shmr.set_pic(ndf, 1, 2, 200.0)

# Close the database
shmr.close_ndf(ndf)