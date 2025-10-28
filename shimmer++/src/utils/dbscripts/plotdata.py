import shimmer as shmr
import matplotlib.pyplot as plt

ndfname = "refined_gasco_autogen.db"

ndf = shmr.open_ndf(ndfname)
cur = ndf.cursor()

nodes = [3] + list(range(43,25,-1)) + [1]

for node in nodes:
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

    plt.plot(data[1])

plt.show()

shmr.close_ndf(ndf)
