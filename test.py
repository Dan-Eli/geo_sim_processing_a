from shapely.geometry import LineString, Point
from shapely.strtree import STRtree
#import random, time
#from rtree import index
from cmath import rect, phase
from math import radians, degrees

def mean_angle(deg):
    a = None
#    a = sum(rect(1, radians(d)) for d,l in deg)
    for d,ll in deg:
        if a is None:
            a = rect(ll, radians(d))
        else:
            a +=  rect(ll, radians(d))

    b = phase(a)
    c = degrees(b)

    d =  degrees(phase(sum(rect(l, radians(d)) for d,l in deg)/len(deg)))

    return (c,d)

def mean_angle2(degrees):

    angle_sum = 0.
    tot_len = 0

    a = sum(rect(l, radians(d)) for d,l in degrees)
    a_phase = phase(a)
    a_degree = degrees(a_phase)


    for deg, len in degrees:
        angle_sum += rect(1, radians(deg))
        tot_len += len

    average_sum = degrees(phase(angle_sum))
    average_sum = average_sum / tot_len

    d = degrees(angle_sum)

    return d



    return degrees(phase(sum(rect(1, radians(d)*l) for d,l in deg)/sum([c[1] for c in deg])))

for angles in [[(350,1000), (10,1)], [(90,1), (180,1), (270,1), (360,1)], [(10,10), (20,1), (30,1)]]:
    print('The mean angle of', angles, 'is:', round(mean_angle(angles)[0], 12), 'degrees')

#for angles in [[(350,2), (10,4)], [(90,2), (180,2), (270,2), (360,2)], [(10,1), (20,2), (30,3)]]:
#    print('The mean angle of', angles, 'is:', round(mean_angle2(angles), 12), 'degrees')

0/0

for xy in [(1,.1),(1,1),(0.1,1),(-0.1,1),(-1,1),(-1,-1),(1,-1)]:
    line0 = LineString([(0,0), xy])
    line1 = LineString([xy, (0,0)])
    for line in (line0,line1):
        x0, y0 = line.coords[0][0], line.coords[0][1]
        x1, y1 = line.coords[1][0], line.coords[1][1]
        delta_y = (y1 - y0)
        delta_x = (x1 - x0)
        angle = math.atan(delta_y / delta_x)
        angle = math.degrees(angle)
        print (x0, y0, x1, y1, angle)

0/0

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
