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

def generate_geojson():
    conn = sqlite3.connect("/home/geoscore/Desktop/GEO++/shimmer/matlab/sqlite/graphs/test_inrete/test_inrete.db")
    cur = conn.cursor()
    cur.execute("SELECT s_number, s_name, s_latitude, s_longitude FROM stations")

    features = []
    for row in cur.fetchall():
        easting, northing = row[2], row[3]  # UTM coordinates
        lon, lat = transformer.transform(easting, northing)  # Convert to WGS84
        features.append({
            "type": "Feature",
            "properties": {"id": row[0], "name": row[1]},
            "geometry": {"type": "Point", "coordinates": [lon, lat]}
        })

    conn.close()

    geojson_data = {"type": "FeatureCollection", "features": features}

    with open("./geojson_data.json", "w") as f:
        json.dump(geojson_data, f)

    print("GeoJSON file generated: geojson_data.json")

if __name__ == "__main__":
    generate_geojson()

