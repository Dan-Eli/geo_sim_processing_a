from shapely.geometry import Point, LineString
from shapely.strtree import STRtree

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

