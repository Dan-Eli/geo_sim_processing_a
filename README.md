# geo_sim_processing

geo_sim_processing is a QGIS plugin that aims to simplify/generalize line and polygon features. It is composed of 3 tools: [Reduce Bend](#Reduce-Bend), [Chordal Axis](#Chordal-Axis) and [Simplifier](#Simplifier)

## Requirements  
- [QGIS](www.qgis.org) >3.14

## QGIS plugin installation
From the GitHub repo download the zip file of the latest tag (or the tag you whish to install) and unzip the content in the QGIS plugin directory _geo_sim_processing_ and reload the plugin geo_sim_processing.  If the _Plugin Reloader_ is not present load it from the menu Pulgins > Manage and Install Plugins 

Plugin directory in Linux: /home/_usename_/.local/share/QGIS/QGIS3/profiles/default/plugins/geo_sim_processing

Plugin directory in Windows: C:\Users\_usename_\AppData\Roaming\QGIS\QGIS3\profiles\default\plugins\geo_sim_processing

Note: Other locations are possible but these are the default one

# Reduce Bend

Reduce Bend is a geospatial simplification and generalization tool for lines and polygons.  Reduce Bend is an implementation and an improvement of the algorithm described in the paper "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this algorithm is that for each line it analyzes its bends (curves) and decides which one to simplify, trying to emulate what a cartographer would do manually to simplify or generalize a line.  Reduce Bend will accept lines and polygons as input.

## Usage

Reduce Bend is a processing script dicoverable in the QGIS Processing Tool Box under Geo Simplification

**Input vector layer**:   Input vector feature to simplify (LineString or Polygon)

**Diameter tolerance**:   Diameter of the minimum adjusted area bend to simplify (to remove) in ground units

**Verbose**:              Flag for extra information

**Exclude polygon**:      Exclude (delete) polygons exteriors below the minimum adjusted area (delete also any interior holes if present)

**Exclude hole**:         Exclude (delete) polygon rings (interior holes) below the minimum adjusted area

**Output vector layer**:  Output vector feature simplified

## Line Simplification versus Line Generalization

*Line Simplification* is the process of removing vertices in a line while trying to keep the maximum number of details within the line whereas *Line Generalization* is the process of removing meaningless (unwanted) details in a line usually for scaling down.  The well known Douglas-Peucker algorithm is a very good example of line simplification tool and Reduce Bend falls more in the category of line generalization tools. Keep in mind thay both algorithms can be complementary because Reduce Bemd will not remove unnecessary vertices in the case of very high densities of vertices.  It may be a good idea to use Douglass Peucker before Reduce Bend in the case of very densed geometries.

## How it works

Reduce Bend will simplify (generalize) lines as well as polygons.  Reduce Bend consists of three main steps: detect bends, determine which bends to simplify and preserve the topological (spatial) relationships.  These 3 steps are detailed below.

* __Detecting bends__ -
For each line and ring composing polygon features, Sherbend will detect the position of each bend.  Wang and Müller defined a bend as being the part of a line which contains a number of subsequent vertices with the inflection angles on all vertices being in opposite sign.
Figure 1a shows a line.  Figure 1b depicts the same line with inflexion signs on ech vertice.  Figure 1c shows the position of the 3 bends each forming an area.

* __Determining the bends to simplify__ -
For each bend of a line or polygon ring, Reduce Bend calculates an adjusted area value using the following formula: *\.75\*A/cmpi* where *A* is the area of the bend *(1)* and *cmpi* the compactness index of the bend.  The compactness index is computed using the following formula: *4\*π\*A/p\*\*2* where *A* is the area of the bend and *p* is the perimeter of the bend. The compactness index varies between \[0..1].  The more circular the bend, the closer the index to 1.  Conversely, the flatter the bend, the closer the index to 0.  The Reduce Bend Diameter tolerance: 4 represents the diameter of a theoretical circle to define the minimum adjusted area value using *\.75\*2\*π\*r\*\*2/cmpi* where *r* is d/2.  Finally, each bend of a line that is below the minimum adjusted area value is replaced by a straight line.  Figure 1d shows the result with the middle bend of the line removed (simplified).

*(1)* The computations are always done in map unit: meters, feet, degrees...

![Figure1](/image/figure1.png)

* __Preserving topological relationship__ -
Before any bend simplifcation is applied, Reduce Bend will always analyze the following 3 topological relationships to ensure they are not affected by the simplification operation: simplicity, intersection and sidedness.  If simplification alters any of those relationships, then it is not performed.  Thereby Sherbend preserves the existing relative topology between the geospatial features to simplify.  

### Simplicity
Reduce Bend will not simplify a bend, if the simplified bend (dashed line in figure 2a) creates a self intersection.  

### Intersection
Reduce Bend will not simplify a bend, if the simplified bend creates an intersection between 2 existing features (figure 2b).  Conflicting features can be a line with another line or a line with a polygon ring.

### Sidedness
Reduce Bend will not simplify a bend, if simplifying the bend creates a sidedness or relative position error between 2 features. Two examples of sidedness issues are shown in figures 2c and 2d.  The preservation of the sidedness topological relationship is particulary important when it comes to simplifying polygon rings.  In figure 2c, simplifying the polygon as shown (dashed line) would make what was an inner hole "pop out" and become external to the new polygon.  In figure 2d, simplifying the bend in the line segment (dashed line) would result in the point feature changing its location relative to the original line.  If for example the original line represents a river and the point represents a building, it would mean the building would find itself on the other side of the river after simplification.   Conflicting features can be a line with a point or a line with a line or a line with a polygon ring.

Note: For any given line or polygon ring, only those bends the simplification of which do not cause any topological issues as expressed above will be simplified.

![Figure2](/image/figure2.png)

### Rule of thumb for the diameter
Reduce Bend can be used for line simplifying often in the context of line generalization. The big question will often be what diameter should we use?  A good starting point is the cartographic rule of thumb -- the *.5mm on the map* -- which says that the minimumm distance between two lines should be greater than 0.5mm on a paper map. So to simplify (generalize) a line for representation at a scale of 1:50 000 for example a diameter of 25m should be a good starting point...

# Chordal Axis

ChordalAxis is a geospatial tool that takes triangles, usually the result of a constraint Delauny trianglulation and creates a skeleton (the center line).  ChordalAxis is an improvement of the algorithm based of the paper "Rectification of the Chordal Axis Transform and a New Criterion for Shape
Decomposition", Lakshman Prasad, 2005".

## Medial Axis Versus Chordal Axis

The skeleton (center line) is a linear feature representation of a polygonized feature. In computational geometry, it is known as the medial axis and many algorithms are approximating it very well.  A major issue with those algorithms is the possible instability for very irregular complex polygons such as dense river or road network polygons. (Figure 4).  The Chordal Axis has shown excellent stability in very very polygons while extracting a very representative skeleton.

## Usage

Chordal Axis is a processing script dicoverable in the QGIS Processing Tool Box under Geo Simplification

**Input vector layer**:   Input vector feature to create the chordal axis (skeleton)

**Correction**:           Flag to correct the skeleton for small centre line, T junction and X junction. Usefull in the case of long any narrow polygon. 

**Output vector layer**:  Output vector feature for the skeleton

 ## How it works

A user will probably create the triangulation from a set of polygons using a constraints Delaunay triangulation tool.  Delaunay triangulation is known to describe polygons well and to  be very robust and stable.  The resulting triangles are the input for the Chordal Axis program.  The Chordal Axis alogorithm will analyze each triangle, determine its type based on the number of adjacent triangles and build the appropriate skeleton (centre line).  All triangles fall within one of the following four types: 1)  _isolated triangle_, when a triangle has no adjacent triangle; 2) _terminal triangle_, when a triangle has only one adjacent triangle; 3) _sleeve triangle_, when a triangle has 2 adjacent triangles; 4) _junction triangle_, when a triangle has 3 adjacent triangles.  Each of the four triangle types will produce a specific centre line.  For the _isolated triangle_, (Figure 3a) no center line (degenerated case) is created; for the _terminal triangle_ (Figure 3b) the mid point of the adjecent side is connected with the opposite angle; for _sleeve triangle_ (Figure 3c) the mid point of the two adajcent sides are connected; for the _junction triangle_ (Figure 3d) the mid points of each side are connected to the centre point of the triangle.  After centre line creation all the centre lines are merged together.  The Chordal Axis transform will preserve [Simplicity](#Simplicity) and [Intersection](#Intersection) topological relationships between the lines forming the skeleton and the outer and inner boundaries of the polygon defined by the Delaynay triangulation.

![figure3](/image/figure3.png)

## Correction
The Chordal Axis algorithm gives a very good approximation of the true medial axis of a polygon but it produces unwanted artifacts when it creates the skeleton especially in the case of long and narrow polygons (figure 5) . The main artifact types are: meaningless small centre line (figure 4a); wrongly formed "T junctions" (figure 4c) and "X junctions" (crossing junction) (figure 4e).  When the correction parameter  is used (-c or --correction), the skeleton will be pruned of the meaningless small centre line (4b); it will correct "T junctions" and rectify the normal direction of the line (figure 4d); and, it will rectify the "X crossing" by merging two T junctions that are adjacent (figure 4f).

![figure5a](/image/figure5a.png "Figure 5a") ![figure5b](/image/figure5b.png "Figure 5b") ![figure5c](/image/figure5c.png "Figure 5c") ![figure5d](/image/figure5d.png "Figure 5d") ![figure5e](/image/figure5e.png "Figure 5e") ![figure5f](/image/figure5f.png "Figure 5f")

  Figure 4a    Figure 4b     Figure 4c    Figure 4d    Figure 4e    Figure 4f

## Rule of thumb for the use of Chordal Axis
Chordal Axis can be used for skeleton extraction and polygon to line transformation in the context of polygon generalization. Often the quality of the skeleton produced will depend on the density of polygon vertices and therefore overall quantity of triangles ingested by Chordal Axis : the more vertices, the higher the number of generated triangles and the better the skeleton (at the price of increased computation time).  Equilateral triangles produce the best skeleton while highly obtuse and/or acute triangles will produce a jagged line that can then be simplified.  The vertex density should not result in either over- or under-simplified features.  Delaunay triangulation and Chordal Axis will give excellent results in very complex situations like a densely polygonized road network such as the one shown in Figure 4 in which all road segments belong to the same polygon!

![figure4](/image/figure4.png)

Figure 5


# Simplifier

Simplify is a geospatial simplification tool for lines and polygons. Simplify implements QGIS's QgsTopologyPreservingSimplifier tool. For line and polygon simplification that tool implements an algorithm similar to the Douglas Peucker algorithm. The implementation preserves the topology within one vector feature but not between vector features. There is also a known bug where the algorithm may create invalid topologies if there are components which are small relative to the tolerance value. In particular, if a small interior hole is very close to an edge, the resulting simplification may result in the hole being moved outside the polygon. This algorithm will detect these situations where one or more rings (interior parts) fall outside the polygon after being simplified and make the polygon invalid. The algoritm will remove (delete) these ring(s) so the feature remains valid after simplification.


Note: While most GIS tools will handle and display invalid geometries like figure 6b, some spatial operation will not be allowed and this is why it's important to keep validity of the geometry after a spatial operation.

![figure6a](/image/figure6a.png "Figure 6a") ![figure6b](/image/figure6b.png "Figure 6b")

        Figure 6a                    Figure 6b


## Usage

Simplifier is a processing script dicoverable in the QGIS Processing Tool Box under Geo Simplification

**Input layer**: The Line String or Polygon layer to simplify

**Tolerance**: The tolerance in ground unit used to simplify the line

**Simplified**: The simplified Line String or Polygon Layer

## How it works

Toposim is an excellent tool to remove vertices on features with high vertex densities.  Try it with small tolerance value and then use [Reduce Bend](#Reduce-Bend) to [generalize features](##Line-Simplification-versus-Line-Generalization).
