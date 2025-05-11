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

def convert_nodes(db_path, features):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT s_number, s_name, s_latitude, s_longitude, s_height FROM stations")

    for row in cur.fetchall():
        lon, lat = transformer.transform(row[2], row[3])
        features.append({
            "type": "Feature",
            "properties": {"id": row[0], "name": row[1], "height": row[4]},
            "geometry": {"type": "Point", "coordinates": [lon, lat]}
        })

    conn.close()

def convert_pipes(db_path, features):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT ROW_NUMBER() OVER() AS NoId, P.p_name, S1.s_latitude, S1.s_longitude, S2.s_latitude, S2.s_longitude FROM pipelines AS P LEFT JOIN stations as S1 ON P.s_from = S1.s_number LEFT JOIN stations as S2 ON P.s_to = S2.s_number")

    for row in cur.fetchall():
        origin_lon, origin_lat = transformer.transform(row[2], row[3])
        dest_lon, dest_lat = transformer.transform(row[4], row[5])
        features.append({
            "type": "Feature",
            "properties": {"id": row[0], "name": row[1]},
            "geometry": {"type": "LineString", "coordinates": [[origin_lon, origin_lat],[dest_lon, dest_lat]]}
        })

    conn.close()

def write_geojson(json_path, features):
    geojson_data = {"type": "FeatureCollection", "features": features}

    with open(json_path, "w") as f:
        json.dump(geojson_data, f)
    print(json_path + " file generated")

def generate_geojson(db_path, json_folder_path):
    nodes = []
    convert_nodes(db_path, nodes)
    write_geojson(json_folder_path + "/nodes.json", nodes)

    pipes = []
    convert_pipes(db_path, pipes)
    write_geojson(json_folder_path + "/pipes.json", pipes)

if __name__ == "__main__":
    db_path = "../graphs/test_inrete/test_inrete.db"
    db_path = "../graphs/test_sicilia/test_sicilia.db"
    json_path = "."
    generate_geojson(db_path, json_path)

