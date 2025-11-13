import shimmer as shmr
import matplotlib.pyplot as plt
import matplotlib.animation as ani

path = ""
ndfname = path + "refined_gasco_autogen.db"

ndf = shmr.open_ndf(ndfname)
cur = ndf.cursor()

nodes13 = [3] + list(range(83,45,-1)) + [1]
nodes23 = [3] + list(range(161,123,-1)) + [2]
nodes24 = [4] + list(range(122,84,-1)) + [2]
nodes34 = [3] + list(range(6,44,1)) + [4]
nodes35 = [3] + list(range(162,200,1)) + [5]

cur.execute("select MAX(timestep) from solution_station_molarfrac");
timesteps = cur.fetchone()[0]

fig, ax = plt.subplots()
line13, = ax.plot( [0]*len(nodes13) )
line23, = ax.plot( [0]*len(nodes23) )
line34, = ax.plot( [0]*len(nodes34) )
line35, = ax.plot( [0]*len(nodes35) )
line24, = ax.plot( [0]*len(nodes24) )
ax.set(ylim = (-0.2,1.2))

values13 = [0] * len(nodes13)
values23 = [0] * len(nodes23)
values24 = [0] * len(nodes24)
values34 = [0] * len(nodes34)
values35 = [0] * len(nodes35)

def animate(timestep):
    print(timestep, timesteps)
    cur.execute(
        """
        select s_number, molarfrac from solution_station_molarfrac
            where timestep = ? and g_name = 14
            order by s_number
        """, (timestep,)
    )
    data = cur.fetchall()

    for i in range(len(nodes13)):
        node = nodes13[i];
        ndata = data[node-1]
        values13[i] = ndata[1]
    line13.set_ydata(values13);

    for i in range(len(nodes23)):
        node = nodes23[i];
        ndata = data[node-1]
        values23[i] = ndata[1]
    line23.set_ydata(values23);

    for i in range(len(nodes24)):
        node = nodes24[i];
        ndata = data[node-1]
        values24[i] = ndata[1]
    line24.set_ydata(values24);

    for i in range(len(nodes34)):
        node = nodes34[i];
        ndata = data[node-1]
        values34[i] = ndata[1]
    line34.set_ydata(values34);

    for i in range(len(nodes35)):
        node = nodes35[i];
        ndata = data[node-1]
        values35[i] = ndata[1]
    line35.set_ydata(values35);

    return line13,line23,line34,line35,line24,

anim = ani.FuncAnimation(fig, animate, interval=1, blit=True,save_count=timesteps)
anim.save("gas.mp4")
#plt.show();


shmr.close_ndf(ndf)
