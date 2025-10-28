import shimmer as shmr
import graphviz

ndfname = "refined.db"


ndf = shmr.open_ndf(ndfname)

dot = graphviz.Graph()

cur = ndf.cursor()
for row in cur.execute("select * from stations"):
    dot.node(str(row[0]))

for row in cur.execute("select * from pipelines"):
    dot.edge(str(row[1]), str(row[2]))


shmr.close_ndf(ndf)

print(dot)