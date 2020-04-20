"""Simplest case to test the Chordal axis"""

from shapely.geometry import LineString

a = LineString([(0,0), (1,1), (2,2), (0,0)])
b = LineString([(10,10), (11,11), (12,12), (10,10)])
case0 = [a,b]

a = LineString([(0,0), (1,1), (2,0), (0,0)])
b = LineString([(1,1), (3,1), (2,0), (1,1)])
c = LineString([(10,10), (11,11), (12,10), (10,10)])
d = LineString([(11,11), (13,11), (12,10), (11,11)])
case1 = [a,b,c,d]

a = LineString([(0,0), (1,1), (2,0), (0,0)])
b = LineString([(1,1), (3,1), (2,0), (1,1)])
c = LineString([(2,0), (3,1), (4,0), (2,0)])

d = LineString([(10,10), (11,11), (12,10), (10,10)])
e = LineString([(11,11), (13,11), (12,10), (11,11)])
f = LineString([(12,10), (13,11), (14,10), (12,10)])
case2 = [a,b,c,d,e,f]

a = LineString([(0,0), (1,1), (2,0), (0,0)])
b = LineString([(1,1), (3,1), (2,0), (1,1)])
c = LineString([(2,0), (3,1), (4,0), (2,0)])
d = LineString([(1,1), (2,2), (3,1), (1,1)])

e = LineString([(10,10), (11,11), (12,10), (10,10)])
f = LineString([(11,11), (13,11), (12,10), (11,11)])
g = LineString([(12,10), (13,11), (14,10), (12,10)])
h = LineString([(11,11), (12,12), (13,11), (11,11)])

case3 = [a,b,c,d, e,f,g,h]

#lst_triangles = case3


ca = ChordalAxis(lst_triangles, GenUtil.ZERO)
#if command.correct:
ca.correct_skeleton()
centre_lines = ca.get_skeleton()