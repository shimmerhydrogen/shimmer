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


###########################################################################
# Shimmer++ tutorial - Programmatically fill a NDF with a sample network
#
# This tutorial will show a possible way to fill the Shimmer++ Network
# Data Files via a Python script.

# Import the 'shimmer' module implemented in 'shimmer.py': the role of
# shimmer.py is just to wrap the SQL statements that will fill the
# tables in the SQL DB. See 'shimmer.py' for details
import shimmer as shmr

# Set a variable with the name of the NDF to create
ndfname = "gasco_autogen.db"

# Create a new NDF
shmr.create_ndf(ndfname)

# Open the NDF just created
ndf = shmr.open_ndf(ndfname)

# In the example network there are two entry stations and two exit stations.
# We want to set some operational limits for those stations; we store them
# in dictionaries.

# A set of limits that we will apply to ENTRY stations
entry_limits = {}
entry_limits['Lmax'] = -300.0
entry_limits['Lmin'] = -10.0
entry_limits['Pmin'] = 6000000.0
entry_limits['Pmax'] = 8000000.0

# A set of limits that we will apply to EXIT stations
exit_limits = {}
exit_limits['Lmin'] = 10.0
exit_limits['Lmax'] = 300.0
exit_limits['Pmin'] = 6000000.0
exit_limits['Pmax'] = 8000000.0

# ------------------------------------------------------------------
# Here we start defining the stations: ReMi stations inject gas in
# the network, and they can inject with constant pressure or with
# a time-varying pressure profile. The profile is specified as a
# list of (time, pressure) pairs. If the profile has only one entry
# the pressure is held constant.

# Station 1
remi_profile_1 = [
    (0, 15914750.0)
]

# Add the station. Parameters:
# NDF file handle   -> ndf
# station number    -> 1
# station type      -> shmr.STATION_REMI
# station height    -> 0.0
shmr.add_station(ndf, 1, shmr.STATION_REMI, 0.0)

# Specify the operational limits of the station. Parameters:
# NDF file handle   -> ndf
# station number    -> 1
# limits dictionary -> entry_limits
shmr.add_remi_limit(ndf, 1, entry_limits)

# Specify the pressure profile for the station. Parameters:
# NDF file handle   -> ndf
# station number    -> 1
# profile           -> remi_profile_1
shmr.add_remi_profile(ndf, 1, remi_profile_1)

# We do exactly the same for Station 2
remi_profile_2 = [
    (0, 16700000.0)
]
shmr.add_station(ndf, 2, shmr.STATION_REMI, 0.0)
shmr.add_remi_limit(ndf, 2, entry_limits)
shmr.add_remi_profile(ndf, 2, remi_profile_2)

# Station 3 is a junction, the only thing to specify is height
shmr.add_station(ndf, 3, shmr.STATION_JUNC, 0.0)

# Station 4
# This is a consumption profile, pairs are (time, flow rate)
cons_profile_4 = [
    (0, 480.0)
]
# Same logic as above: add station 4 of shmr.STATION_CONS type
shmr.add_station(ndf, 4, shmr.STATION_CONS, 0.0)
# Set operational limits
shmr.add_cons_limit(ndf, 4, exit_limits)
# Set mass flow profile
shmr.add_cons_profile(ndf, 4, cons_profile_4)

# Same as above for Station 5
cons_profile_5 = [
    (0, 72.0)
]
shmr.add_station(ndf, 5, shmr.STATION_CONS, 0.0)
shmr.add_cons_limit(ndf, 5, exit_limits)
shmr.add_cons_profile(ndf, 5, cons_profile_5)

# ------------------------------------------------------------------
# Now the pipes. We are going to do quality tracking, so we need to specify
# the number of segments in which the pipes are discretized
refinement_segs = 40

# Add to the ndf a pipe from node 3 to node 4
shmr.add_pipe(ndf, 3, 4)
# For the pipe above, set the diameter (0.9644 m), the length (619000.0 m),
# the roughness (1.4e-5) and the number of segments for the discretization
shmr.set_pipe_params(ndf, 3, 4, 0.9664, 619000.0, 1.4e-05, refinement_segs)

# Same as above for the other pipes
shmr.add_pipe(ndf, 1, 3)
shmr.set_pipe_params(ndf, 1, 3, 0.9664, 300000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 2, 4)
shmr.set_pipe_params(ndf, 2, 4, 1.016, 658000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 2, 3)
shmr.set_pipe_params(ndf, 2, 3, 0.6602, 228000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 3, 5)
shmr.set_pipe_params(ndf, 3, 5, 1.016, 840000.0, 1.4e-05, refinement_segs)


# ------------------------------------------------------------------
# To start a simulation we need to specify the initial condition of
# each station. Here we create a list of triples with
# (station number, pressure initial condition, flow rate initial condition)
# Initial conditions for stations
station_ics = [
    (1, 15914750.0, -282.4514),
    (2, 16700000.0, -282.4514),
    (3, 15449982.5100397, 1.0e-08),
    (4, 13840652.4843216, 479.667250138719),
    (5, 15135717.6888121, 79.9555703310221)
]

# With a for cycle we iterate on the list above and we call shmr.set_sic()
# to actually set the initial condition
for sic in station_ics:
    shmr.set_sic(ndf, sic[0], sic[1], sic[2])

# Same thing for pipes: the triples in the list are
# (from station, to station, flow rate)
# Initial conditions for pipes
pipe_ics = [
    (3, 4, 182.079365649059),
    (1, 3, 146.915345084878),
    (2, 4, 297.587884489659),
    (2, 3, 115.119590905204),
    (3, 5, 79.9555703310221)
]

# And again we iterate on the list and set the pipe initial conditions
# by calling shmr.set_pic()
for pic in pipe_ics:
    shmr.set_pic(ndf, pic[0], pic[1], pic[2])

# ------------------------------------------------------------------
# Now the gas molar fractions:
# Gas molar fractions: this example only supports 4 gases, it is
# trivial to extend 'shimmer.py' to support all of them

# Dictionary of the molar fractions
gmf1 = {}
gmf1['CH4'] = 0.8   # 80% CH4
gmf1['N2'] = 0.0    #  0% N2
gmf1['CO2'] = 0.0   #  0% CO2
gmf1['H2'] = 0.2    # 20% H2
# assign that molar fraction to station 1
shmr.set_gmf(ndf, 1, gmf1) 

# Same as above for the other stations
gmf2 = {}
gmf2['CH4'] = 0.9
gmf2['N2'] = 0.0
gmf2['CO2'] = 0.0
gmf2['H2'] = 0.1
shmr.set_gmf(ndf, 2, gmf2)

gmf3 = {}
gmf3['CH4'] = 1.0
gmf3['N2'] = 0.0
gmf3['CO2'] = 0.0
gmf3['H2'] = 0.0
shmr.set_gmf(ndf, 3, gmf3)

gmf4 = {}
gmf4['CH4'] = 1.0
gmf4['N2'] = 0.0
gmf4['CO2'] = 0.0
gmf4['H2'] = 0.0
shmr.set_gmf(ndf, 4, gmf4)

gmf5 = {}
gmf5['CH4'] = 1.0
gmf5['N2'] = 0.0
gmf5['CO2'] = 0.0
gmf5['H2'] = 0.0
shmr.set_gmf(ndf, 5, gmf5)
