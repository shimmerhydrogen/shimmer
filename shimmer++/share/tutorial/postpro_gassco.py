import shimmer as shmr
import matplotlib.pyplot as plt
import matplotlib.animation as ani

path = ""
ndfname = path + "refined_gassco_autogen.db"

ndf = shmr.open_ndf(ndfname)
cur = ndf.cursor()

nodes13 = [3] + list(range(83,45,-1)) + [1]


cur.execute("select MAX(timestep) from solution_station_molarfrac");
#cur.execute("select MAX(timestep) from solution_station_pressures");
timesteps = cur.fetchone()[0]

fig, ax = plt.subplots()
line13, = ax.plot( [0]*len(nodes13) )
ax.set(ylim = (-0.2,1.2))
#ax.set(ylim = (0,20000000))

values13 = [0] * len(nodes13)

def animate(timestep):
    print(timestep, timesteps)
    cur.execute(
        """
        select s_number, molarfrac from solution_station_molarfrac
            where g_name = 14 and timestep = ?
            order by s_number
        """, (timestep,)
    )
    
    #cur.execute(
    #    """
    #    select s_number, pressure from solution_station_pressures
    #        where timestep = ?
    #        order by s_number
    #    """, (timestep,)
    #)
    
    data = cur.fetchall()

    for i in range(len(nodes13)):
        node = nodes13[i];
        ndata = data[node-1]
        values13[i] = ndata[1]
    line13.set_ydata(values13);

    return line13,

anim = ani.FuncAnimation(fig, animate, interval=20, blit=True,save_count=timesteps)
anim.save("gas.mp4")
#plt.show();


shmr.close_ndf(ndf)
