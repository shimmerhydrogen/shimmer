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

import os
import sqlite3
def create_ndf(filename):
    if os.path.exists(filename):
        print(f"NDF already exists: {filename}")
        return False

    try:
        with open("shimmer.sql", "r") as file:
            sql_content = file.read()
    except FileNotFoundError:
        print("Could not find database schema")
        return False

    try:
        with sqlite3.connect(filename) as conn:
            cursor = conn.cursor()
            cursor.executescript(sql_content)
            conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
        return False

    return True

STATION_REMI = 1
STATION_INJ = 2
STATION_CONS = 3
STATION_JUNC = 4

def open_ndf(filename):
    try:
        conn = sqlite3.connect(filename)
        return conn
    except sqlite3.Error as e:
        print(f"Error opening database: {e}")
        return None

def close_ndf(conn):
    if conn:
        conn.close()

def add_station(conn, station, type, height):
    try:
        cursor = conn.cursor()

        cursor.execute(
            """
            INSERT INTO stations (s_number, s_name, t_type, s_height)
            VALUES (?, ?, ?, ?)
            """,
            (station, f"station_{station}", type, height)
        )

        #cursor.execute(
        #    """
        #    INSERT INTO gas_molar_fractions (s_number, frac_CH4)
        #    VALUES (?, ?)
        #    """,
        #    (station, 1.0)
        #)

        conn.commit()

    except sqlite3.Error as e:
        print(f"SQLite error in add_station: {e}")
        conn.rollback()

def add_remi_limit(conn, station, limit):
    try:
        cursor = conn.cursor()

        cursor.execute(
            """
            INSERT INTO limits_remi_wo (s_number, lim_Lmin, lim_Lmax, lim_Pmin, lim_Pmax)
            VALUES (?, ?, ?, ?, ?)
            """,
            (station, limit['Lmin'], limit['Lmax'], limit['Pmin'], limit['Pmax'])
        )

        conn.commit()

    except sqlite3.Error as e:
        print(f"SQLite error in add_remi_limit: {e}")
        conn.rollback()

def add_remi_profile(conn, station, profile):
    try:
        cursor = conn.cursor()
        cursor.execute("BEGIN")

        for row in profile:
            cursor.execute(
                "INSERT INTO profiles_remi_wo VALUES (?, ?, ?)",
                (station, row[0], row[1])
            )

        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_remi_profile: {e}")
        conn.rollback()


def add_inj_limit(conn, station, limit):
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO limits_injection_w(s_number, lim_Lmin, lim_Lmax, lim_Pmin, lim_Pmax)
            VALUES (?, ?, ?, ?, ?)
            """,
            (station, limit['Lmin'], limit['Lmax'], limit['Pmin'], limit['Pmax'])
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_inj_limit: {e}")
        conn.rollback()


def add_inj_profile(conn, station, profile):
    try:
        cursor = conn.cursor()
        cursor.execute("BEGIN")

        for row in profile:
            cursor.execute(
                "INSERT INTO profiles_injection_w VALUES (?, ?, ?, ?)",
                (station, row[0], row[1], row[2])
            )

        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_inj_profile: {e}")
        conn.rollback()


def add_cons_limit(conn, station, limit):
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO limits_conspoint_wo VALUES (?, ?, ?, ?, ?)",
            (station, limit['Lmin'], limit['Lmax'], limit['Pmin'], limit['Pmax'])
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_cons_limit: {e}")
        conn.rollback()


def add_cons_profile(conn, station, profile):
    try:
        cursor = conn.cursor()
        cursor.execute("BEGIN")

        for row in profile:
            cursor.execute(
                "INSERT INTO profiles_conspoint_wo VALUES (?, ?, ?)",
                (station, row[0], row[1])
            )

        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_cons_profile: {e}")
        conn.rollback()


def set_sic(conn, station, P, L):
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO station_initial_conditions VALUES (?, ?, ?)",
            (station, P, L)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_sic: {e}")
        conn.rollback()


def add_pipe(conn, from_station, to_station):
    pname = f"pipe_{from_station}_{to_station}"
    PLAIN_PIPE = 0
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO pipelines(p_name, s_from, s_to, p_type) VALUES (?, ?, ?, ?)",
            (pname, from_station, to_station, PLAIN_PIPE)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_pipe: {e}")
        conn.rollback()


def set_pipe_params(conn, from_station, to_station, diam, length, rough, ref_nsegs = 0):
    pname = f"pipe_{from_station}_{to_station}"
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO pipe_parameters VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (pname, from_station, to_station, diam, length, rough, ref_nsegs)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_pipe_params: {e}")
        conn.rollback()

def add_compressor(conn, from_station, to_station):
    pname = f"pipe_{from_station}_{to_station}"
    COMPR_STAT = 1
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO pipelines(p_name, s_from, s_to, p_type) VALUES (?, ?, ?, ?)",
            (pname, from_station, to_station, COMPR_STAT)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in add_compressor: {e}")
        conn.rollback()

ON_POWER    =  0
ON_OPRESS   =  1
ON_IPRESS   =  2
ON_RATIO    =  3
ON_MASSFLOW =  4
OFF_BYPASS  = 10
OFF_CLOSED  = 11

def set_compressor_limits(conn, from_station, to_station, limits):
    pname = f"pipe_{from_station}_{to_station}"
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO compressor_limits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (pname, from_station, to_station, limits['max_power'],
                limits['max_outpress'], limits['min_inpress'],
                limits['max_ratio'], limits['min_ratio'],
                limits['max_massflow']
            )
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_compressor_limits: {e}")
        conn.rollback()

def set_pic(conn, from_station, to_station, G):
    pname = f"pipe_{from_station}_{to_station}"
    try:
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO pipe_initial_conditions VALUES (?, ?, ?, ?)",
            (pname, from_station, to_station, G)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_pic: {e}")
        conn.rollback()

def set_gmf(conn, station, gmf):
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO gas_molar_fractions(s_number,frac_CH4,
                frac_N2,frac_CO2,frac_H2) VALUES (?, ?, ?, ?, ?)
            """,
            (station, gmf['CH4'], gmf['N2'], gmf['CO2'], gmf['H2'])
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_gmf: {e}")
        conn.rollback()
