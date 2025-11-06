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

ndfname = "gasco_autogen.db"

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

# ------------------------------------------------------------------
# Station 1
remi_profile_1 = [
    (0, 15914750.0)
]
shmr.add_station(ndf, 1, shmr.STATION_REMI, 0.0)
shmr.add_remi_limit(ndf, 1, entry_limits)
shmr.add_remi_profile(ndf, 1, remi_profile_1)

# Station 2
remi_profile_2 = [
    (0, 16700000.0)
]
shmr.add_station(ndf, 2, shmr.STATION_REMI, 0.0)
shmr.add_remi_limit(ndf, 2, entry_limits)
shmr.add_remi_profile(ndf, 2, remi_profile_2)

# Station 3
shmr.add_station(ndf, 3, shmr.STATION_JUNC, 0.0)

# Station 4
cons_profile_4 = [
    (0, 480.0)
]
shmr.add_station(ndf, 4, shmr.STATION_CONS, 0.0)
shmr.add_cons_limit(ndf, 4, exit_limits)
shmr.add_cons_profile(ndf, 4, cons_profile_4)

# Station 5
cons_profile_5 = [
    (0, 72.0)
]
shmr.add_station(ndf, 5, shmr.STATION_CONS, 0.0)
shmr.add_cons_limit(ndf, 5, exit_limits)
shmr.add_cons_profile(ndf, 5, cons_profile_5)

# ------------------------------------------------------------------
refinement_segs = 40
shmr.add_pipe(ndf, 3, 4)
shmr.set_pipe_params(ndf, 3, 4, 0.9664, 619000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 1, 3)
shmr.set_pipe_params(ndf, 1, 3, 0.9664, 300000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 2, 4)
shmr.set_pipe_params(ndf, 2, 4, 1.016, 658000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 2, 3)
shmr.set_pipe_params(ndf, 2, 3, 0.6602, 228000.0, 1.4e-05, refinement_segs)

shmr.add_pipe(ndf, 3, 5)
shmr.set_pipe_params(ndf, 3, 5, 1.016, 840000.0, 1.4e-05, refinement_segs)


# ------------------------------------------------------------------
# Initial conditions for stations
station_ics = [
    (1, 15914750.0, -282.4514),
    (2, 16700000.0, -282.4514),
    (3, 15449982.5100397, 1.0e-08),
    (4, 13840652.4843216, 479.667250138719),
    (5, 15135717.6888121, 79.9555703310221)
]

for sic in station_ics:
    shmr.set_sic(ndf, sic[0], sic[1], sic[2])

# Initial conditions for pipes
pipe_ics = [
    (3, 4, 182.079365649059),
    (1, 3, 146.915345084878),
    (2, 4, 297.587884489659),
    (2, 3, 115.119590905204),
    (3, 5, 79.9555703310221)
]

for pic in pipe_ics:
    shmr.set_pic(ndf, pic[0], pic[1], pic[2])

# ------------------------------------------------------------------
# Gas molar fractions
gmf1 = {}
gmf1['CH4'] = 0.8
gmf1['N2'] = 0.0
gmf1['CO2'] = 0.0
gmf1['H2'] = 0.2
shmr.set_gmf(ndf, 1, gmf1)

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
