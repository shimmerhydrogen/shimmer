import os
import sqlite3
def create_ndf(filename):
    if os.path.exists(filename):
        print(f"NDF already exists: {filename}")
        return False

    try:
        with open("../../../../sqlite/shimmer.sql", "r") as file:
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

        cursor.execute(
            """
            INSERT INTO gas_molar_fractions (s_number, frac_CH4)
            VALUES (?, ?)
            """,
            (station, 1.0)
        )

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


def set_pipe_params(conn, from_station, to_station, diam, length, rough):
    pname = f"pipe_{from_station}_{to_station}"
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO pipe_parameters(p_name, s_from, s_to, diameter, length, roughness)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (pname, from_station, to_station, diam, length, rough)
        )
        conn.commit()
    except sqlite3.Error as e:
        print(f"SQLite error in set_pipe_params: {e}")
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
