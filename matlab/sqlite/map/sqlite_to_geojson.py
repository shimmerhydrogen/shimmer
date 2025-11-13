import sqlite3
import json
from pyproj import Transformer

# EPSG Code	Name	Type	Usage
# 4326	WGS 84	Geographic (Lat/Lon)	Global GPS, OpenStreetMap, Web mapping
# 3857	Pseudo-Mercator	Web Mercator (Meters)	Web maps (Google Maps, OpenLayers, Leaflet)
# 25833	ETRS89 / UTM Zone 33N	Projected (Meters)	European standard for mapping
# 3004	Monte Mario / Italy Zone 2	Projected (Meters)	Older Italian cartographic system
# 7794	RDN2008 / UTM Zone 33N	Projected (Meters)	Modern Italian reference system
# 6706	RDN2008	Geographic (Lat/Lon)	Replacement for Monte Mario

transformer = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True) # EU standard
#transformer = Transformer.from_crs("EPSG:4326", "EPSG:4326", always_xy=True) # WGS 84

def chk_conn(db_path):
     try:
        conn = sqlite3.connect(db_path)
        conn.cursor()
        conn.close()
        return True
     except Exception as ex:
        return False

def convert_nodes(db_path, features, nodes_map):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT s_number, s_name, s_latitude, s_longitude, s_height, t_type FROM stations")

    node_index = 0
    for row in cur.fetchall():
        node_id = row[0]
        lon, lat = transformer.transform(row[2], row[3])
        features.append({
            "type": "Feature",
            "properties": {"id": node_id, "name": row[1], "height": row[4], "type": row[5], "pressure": []},
            "geometry": {"type": "Point", "coordinates": [lon, lat]}
        })
        nodes_map[node_id] = node_index
        node_index = node_index + 1

    conn.close()
    print('Imported', node_index, 'nodes')

def convert_nodes_solution(db_path, features, nodes_map, time_step):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT s_number, timestep, pressure FROM solution_station_pressures WHERE timestep = {0}".format(time_step))

    for row in cur.fetchall():
        node_id = row[0]
        node_index = nodes_map[node_id]
        features[node_index]["properties"]["pressure"] = row[2]
    
    cur = conn.cursor()
    cur.execute("SELECT SM.s_number, SM.timestep, G.g_formula, SM.molarfrac FROM solution_station_molarfrac AS SM LEFT JOIN gases as G ON SM.g_name = G.g_num WHERE molarfrac > 1e-4 AND timestep = {0}".format(time_step))    

    for row in cur.fetchall():
        node_id = row[0]
        node_index = nodes_map[node_id]
        gas_number = row[2]
        features[node_index]["properties"]["gas_" + gas_number] = row[3]        

    conn.close()

def convert_pipes(db_path, features, pipes_map):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT ROW_NUMBER() OVER() AS NoId, P.p_name, S1.s_latitude, S1.s_longitude, S2.s_latitude, S2.s_longitude, P.p_type FROM pipelines AS P LEFT JOIN stations as S1 ON P.s_from = S1.s_number LEFT JOIN stations as S2 ON P.s_to = S2.s_number")

    pipe_index = 0
    for row in cur.fetchall():
        pipe_id = row[0]
        origin_lon, origin_lat = transformer.transform(row[2], row[3])
        dest_lon, dest_lat = transformer.transform(row[4], row[5])
        features.append({
            "type": "Feature",
            "properties": {"id": pipe_id, "name": row[1], "type": row[6], "flowrate": []},
            "geometry": {"type": "LineString", "coordinates": [[origin_lon, origin_lat],[dest_lon, dest_lat]]}
        })
        pipes_map[pipe_id] = pipe_index
        pipe_index = pipe_index + 1

    conn.close()
    print('Imported', pipe_index, 'pipes')

def convert_pipes_solution(db_path, features, pipes_map, time_step):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT ROW_NUMBER() OVER() AS NoId, timestep, flowrate FROM solution_pipe_flowrates WHERE timestep = {0}".format(time_step))

    for row in cur.fetchall():
        pipe_id = row[0]
        pipe_index = pipes_map[pipe_id]
        features[pipe_index]["properties"]["flowrate"] = row[2]

    conn.close()

def write_geojson(json_path, features):
    geojson_data = {"type": "FeatureCollection", "features": features}

    with open(json_path, "w") as f:
        json.dump(geojson_data, f)
    print(json_path + " file generated")

def generate_geojson(db_path, json_folder_path, time_step):
    nodes = []
    nodes_map = {}
    convert_nodes(db_path, nodes, nodes_map)
    convert_nodes_solution(db_path, nodes, nodes_map, time_step)
    write_geojson(json_folder_path + "/nodes_T{0}.json".format(time_step), nodes)

    pipes = []
    pipes_map = {}
    convert_pipes(db_path, pipes, pipes_map)
    convert_pipes_solution(db_path, pipes, pipes_map, time_step)
    write_geojson(json_folder_path + "/pipes_T{0}.json".format(time_step), pipes)

if __name__ == "__main__":
    db_path = "./inRete.db"
    json_path = "."

    if not chk_conn(db_path):
        print("ERROR on connection to DB {0}".format(db_path))
        exit()
    
    generate_geojson(db_path, json_path, 0)
    generate_geojson(db_path, json_path, 100)
    generate_geojson(db_path, json_path, 250)
    generate_geojson(db_path, json_path, 500)
    generate_geojson(db_path, json_path, 1000)

