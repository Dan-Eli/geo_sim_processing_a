from shapely.geometry import Point, LineString
from shapely.strtree import STRtree
from lib_geosim import SpatialContainer, LineStringSc
import random


lst_lines = []

for i in range(100000):
    offset = 10000.
    x = random.random() * offset
    y = random.random() * offset
    coords = [(x,y),(x+5, y+5),(x,y+10),(x,y)]
    lst_lines.append(LineStringSc(coords))

s_cont = SpatialContainer()
s_cont.add_features(lst_lines)
for i in range(100000):
    box = s_cont.get_features(Point(5000,5000).bounds)
    box = s_cont.get_features(Point(5000,5000).bounds)
    box = s_cont.get_features(Point(5000,5000).bounds)
    if i%1000 == 0:
        print (i)
0/0


tree = STRtree(lst_lines)
for i in range(100000):
    box = tree.query(Point(5000,5000))
    box = tree.query(Point(5000, 5000))
    box = tree.query(Point(5000, 5000))
    if i%1000 == 0:
        print (i)

print ('fin')
0/0




tree = STRtree([])
box = tree.query(Point(4,4).buffer(0.99))
print (len(box))
lines = []
for i in range(10):
    line = LineString([(0,0),(i,i)])
    line.__att = i+14
    lines.append(line)

tree = STRtree(lines)
line = lines[5]
line.__att1 = "coco"
line.__att = "coco"
box = tree.query(Point(4,4).buffer(0.99))

print (len(box))
for line in box:
    print (line.__att)
    line._att1 = 27


print ('fin')

