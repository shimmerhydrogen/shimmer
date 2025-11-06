import shimmer as shmr
import matplotlib.pyplot as plt

path = ""
ndfname = path + "refined_gasco_autogen.db"

ndf = shmr.open_ndf(ndfname)
cur = ndf.cursor()

nodes = [3] + list(range(43,25,-1)) + [1]

num_lines = len(nodes)
colors  = plt.get_cmap('jet', num_lines) 
markers = ['+','x','o','*','s','^','v','D','>','<','p','h']


for i in range (num_lines):
    node = nodes[i]
    print(node)
    cur.execute(
        """
        select timestep, molarfrac from solution_station_molarfrac
            where s_number = ? and g_name = 2
            order by timestep
        """, (node,)
    )
    data = cur.fetchall()
    data = list(zip(*data))

    clr = colors(i % colors.N)
    mkr = markers[i % len(markers)]
    
    plt.plot(data[1], color=clr,  marker=mkr, markevery=20, label="St."+str(node))

plt.title("Molar fractions")
plt.ylabel("y [-]")        
plt.xlabel("Time [s]")
plt.legend()
plt.show()

shmr.close_ndf(ndf)
