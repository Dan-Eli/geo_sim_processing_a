from shapely.geometry import LineString
from shapely.strtree import STRtree
import random, time
from rtree import index

lst_lines = []
lst_intersects = []

# Create the triangles
for i in range(250000):
    x = random.random() * 10000.
    y = random.random() * 10000.
    coords = [(x,y),(x+5, y+5),(x,y+10),(x,y)]
    lst_lines.append(LineString(coords))

# Create the bounding boxes
for i in range(10000):
    x = random.random() * 10000.
    y = random.random() * 10000.
    coords = [(x,y),(x+15,y),(x+15,y+15),(x,y+15),(x,y)]
    lst_intersects.append(LineString(coords))

# Create shapely STRtree
tree = STRtree(lst_lines)

# Create RTree
idx = index.Index()
for i, line in enumerate(lst_lines):
    idx.insert(i, line.bounds)
print (time.time())

sec1 = time.time()

# finf the intersection with STRtree
str_tree_nbr = 0
for intersect in lst_intersects:
    str_tree = tree.query(intersect)
    str_tree_nbr += len(str_tree)

sec2 = time.time()
print("Seconds for STRtree =", sec2-sec1)
print ("Str tree number: ", str_tree_nbr)

# Find the intersections with RTree
rtree_nbr = 0
for intersect in lst_intersects:
    rtree = idx.intersection(intersect.bounds)
    rtree_nbr += len(list(rtree))

sec3 = time.time()
print("Seconds for RTree =", sec3-sec2)
print ("Rtree number: ", rtree_nbr)
