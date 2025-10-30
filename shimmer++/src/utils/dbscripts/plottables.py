import shimmer as shmr
import matplotlib.pyplot as plt

def find_stations(cur, init_st, output_st):

    cur.execute("SELECT s_number, s_name FROM stations WHERE s_name LIKE ?", (f'station_{init_st}%',))
    stations = cur.fetchall()

    pattern = f"fict_({init_st}_{output_st})_%"
    cur.execute(
        "SELECT s_number, s_name FROM stations WHERE s_name LIKE ?",
        (pattern,)
    )
    stations.extend(cur.fetchall())

    cur.execute("SELECT s_number, s_name FROM stations WHERE s_name LIKE ?", (f'station_{output_st}%',))
    stations.extend(cur.fetchall())

    if not stations:
        print("No stations found for the given init/output points.")
        ndf.close()
        exit()
    
    return stations

def plot_var(stations, var_tab, var_name):

    table_name = f"solution_station_{var_tab}"
    plt.figure(figsize=(8, 5))
    
    colors  = plt.get_cmap('jet', len(stations)) 
    markers = ['+','x','o','*','s','^','v','D','>','<','p','h']
    
    i = 0

    for s_number, s_name in stations:

        cur.execute(
            f"""
            SELECT timestep, {var_name} FROM {table_name} 
                WHERE s_number = ?
                ORDER BY timestep
            """,(s_number,)
        )
        data = cur.fetchall()

        if not data:
            print(f"No data found in {table_name} for station {s_number}")
            continue

        timesteps = [row[0] for row in data]
        values = [row[1] for row in data]

        clr = colors(i % colors.N)
        mkr = markers[i % len(markers)]
       
        plt.plot(timesteps, values, color=clr,  marker=mkr, markevery=20, label=s_name)
        
        i=i+1
    plt.title(f"{var_name.capitalize()} vs Timestep")
    plt.xlabel("Time steps")
    plt.ylabel(var_name.capitalize())
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
def plot_molarfrac(stations, g_number):
    table_name = f"solution_station_molarfrac"
    plt.figure(figsize=(8, 5))
    
    colors  = plt.get_cmap('jet', len(stations)) 
    markers = ['+','x','o','*','s','^','v','D','>','<','p','h']
    
    i=0
    for s_number, s_name in stations:
        cur.execute(
            f"""
            select timestep, molarfrac from solution_station_molarfrac
            where s_number = ? and g_name = ?
            order by timestep
            """,(s_number,g_number)
        )
        data = cur.fetchall()
        data = list(zip(*data))

        if not data:
            print(f"No data found in {table_name} for station {s_number}")
            continue

        timesteps = data[0]
        values = data[1]

        clr = colors(i % colors.N)
        mkr = markers[i % len(markers)]
       
        plt.plot(timesteps, values, color=clr,  marker=mkr, markevery=20, label=s_name)
        
        i=i+1
    plt.title(f"Molar fractions vs Timestep")
    plt.xlabel("Time steps")
    plt.ylabel("Molar fraction")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

#------------------------------------------------------------------------------
#  MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN   MAIN
#------------------------------------------------------------------------------

init_st = 1     
output_st = 3   
path = ""
ndfname = path + "refined_gasco_autogen.db"

# --- CONNECT TO DATABASE ---
ndf = shmr.open_ndf(ndfname)
cur = ndf.cursor()


stations  = find_stations(cur, init_st, output_st)
variables = {"pressures": "pressure", "flowrates":"flowrate"}

plot_molarfrac(stations, 0)
for var_tab, var_name in variables.items():     
    plot_var(stations, var_tab, var_name)

ndf.close()

