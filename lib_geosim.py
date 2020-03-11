#!/usr/bin/env python
# -=- encoding: utf-8 -=-

"""
General classes and utilities needed for the GeoSim.

"""

import math
from shapely.geometry import Point, LineString, LinearRing, Polygon
from shapely.ops import linemerge
from collections.abc import Iterable
from collections import OrderedDict
import fiona

try:
    from rtree import Rtree

    lib_Rtree = True
except:
    # If Rtree is not installed
    lib_Rtree = False
    from shapely.strtree import STRtree


class LineStringSc(LineString):
    """LineString specialization to be included in the SpatialContainer"""

    def __init__(self, coords):
        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        return super().coords

    @coords.setter
    def coords(self, coords):
        # Update the coord attribute in the parent class
        LineString.coords.__set__(self, coords)

        if self._sc_scontainer != None:  # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_scontainer.update_spatial_index(self)


class PointSc(Point):
    """LineString specialization to be included in the SpatialContainer"""

    def __init__(self, coords):
        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        return super().coords

    @coords.setter
    def coords(self, coords):
        Point.coords.__set__(self, coords)

        if self._sc_scontainer is not None:  # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_container.update_bbox(self)


class PolygonSc(Polygon):
    """Polygon specialization to be included in the SpatialContainer"""

    def __init__(self, exterior, interiors=None):
        super().__init__(exterior, interiors)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def exterior(self):
        return super().exterior

    @property
    def interiors(self):
        return super().interiors

    @exterior.setter
    def exterior(self, exterior):
        raise GeoSimException("Cannot update the exterior coordinates of a polygon")

    @interiors.setter
    def interiors(self, interiors):
        raise GeoSimException("Cannot update the interior coordinates of a polygon")


class GenUtil:
    """This class defines a series of generic utility static method"""

    # Define some constants...
    POINT = 'Point'
    LINE_STRING = 'LineString'
    POLYGON = 'Polygon'
    POLYGON_EXTERIOR = 'PolygonExterior'
    POLYGON_INTERIOR = 'PolygonInterior'

    ANTI_CLOCKWISE = 1
    CLOCKWISE = -1
    CROSSING_LINE = 'Crossing line'
    SIDEDNESS = 'Sidedness'
    SIMPLE_LINE = 'Simple line'
    INVALID = 'Invalid'
    ZERO = 0.000001
    RADIAN = 'Radian'
    DEGREE = 'Degree'

    @staticmethod
    def make_iterable(iter_feature):
        """Test if the parameter is iterable; if not make it iterable by creating a tuple of one element

        *Parameters*:
            - iter_feature: Object to test for iterability

        *Returns*:
            - Iterable object

        """

        if not isinstance(iter_feature, Iterable):
            iter_feature = (iter_feature,)

        return iter_feature

    @staticmethod
    def distance(p1, p2):
        """Calculate the euclidean distance between 2 points

        *Parameters*:
            - p1: (x,y) tuple of the first coordinate
            - p2: (x,y) tuple of the second coordinate

        *Returns*:
            - Distance between the 2 points (real)

        """

        return math.sqrt((p2[0] - p1[0]) ** 2.0 + (p2[1] - p1[1]) ** 2.0)

    @staticmethod
    def compute_angle(p1, p2, p3, type=DEGREE):
        """
        Function to calculate angle between two vectors.
        """

        return GenUtil.angle_vector(p1, p2, p3, type)

    @staticmethod
    def angle_vector(p1, p2, p3, type=DEGREE):
        """Calculate the angle formed by the vector p1-p2 and p2-p3

        *Parameters*:
            - p1: (x,y) tuple of coordinates
            - p2: (x,y) tuple of coordinates
            - p3: (x,y) tuple of coordinates
            - type: Angle type DEGREE or ANGLE

        *Returns*:
            - The angle between the vector p1-p2 and p2-p3 (float)

        """

        a = (p2[0] - p1[0], p2[1] - p1[1])
        b = (p2[0] - p3[0], p2[1] - p3[1])
        len_a = (a[0] ** 2. + a[1] ** 2.) ** .5
        len_b = (b[0] ** 2. + b[1] ** 2.) ** .5

        dot_p = a[0] * b[0] + a[1] * b[1]

        # If P1 == P2 or P2 == P3 ===> angle is 180.
        if len_a * len_b != 0.0:
            value = dot_p / (len_a * len_b)
            if value >= 1.0:
                value = 1.0
            if value <= -1.0:
                value = -1.0
        else:
            value = -1.0

        theta = math.acos(value)

        if type == GenUtil.DEGREE:
            theta = math.degrees(theta)

        return theta

    @staticmethod
    def orientation(p0, p1, p2):
        """ Calculate the orientation (clockwise or anticlockwise) of a line formed by 3 vertices using the dot product

        Parameters:
            p0, p1, p2: Three (x,y) coordinates tuple

        Return value
            float the direction of the line an
                0 : Straight line
                1: Counter clockwise angle
                -1 : Clockwise angle

        """

        orient = ((p0[0] - p1[0]) * (p2[1] - p1[1])) - ((p2[0] - p1[0]) * (p0[1] - p1[1]))

        if orient > 0.:
            orient = 1
        elif orient < 0.:
            orient = -1
        else:
            orient = 0

        return orient

    @staticmethod
    def rescale_vector(p1, p2, scale_factor):
        """Rescale the vector defined by the points P1 and P2 by a factor
        of SCALE_FACTOR

        Rescale the vector defined by the points P1 and P2 by a factor
        of SCALE_FACTOR

        *Parameters*:
            - P1: First coordinate of the vector. Tuple (x,y)
            - P2: Last  coordinate vector to rescale (second point)
            - scale_factor: factor to scale the vector (same for x and y)

        *Returns*: *TBA*

        """

        x1 = p1[0]
        y1 = p1[1]
        x2 = p2[0]
        y2 = p2[1]

        vec_x = x2 - x1
        vec_y = y2 - y1

        vec_x = vec_x * scale_factor
        vec_y = vec_y * scale_factor

        x_out = vec_x + x1
        y_out = vec_y + y1

        return (x_out, y_out)

    @staticmethod
    def mid_point(p1, p2):
        """Return a point in the middle of the 2 points

        *Parameters*:
            - p1: x,y tuple for the first point
            - p2: x,y tuple for the second point

        *Returns*:
            - x,y tuple of the milddle point
        """

        x = (p1[0] + p2[0]) / 2.
        y = (p1[1] + p2[1]) / 2.

        return (x, y)

    @staticmethod
    def calculate_compactness_index(area, perimeter):
        """Calculate the compactness index based of the perimeter and area

        Args:
            area (float): Area of the polygon
            perimeter (float): Perimeter of the area

        Return:
            (float): Compactness index

        """

        return 4 * area * math.pi / (perimeter ** 2.0)

    @staticmethod
    def build_bounding_box(tolerance, coord):
        """Create and adjust a bounding box (xmin, ymin, xmax, ymax) with a small tolerance

        *Parameters*:
            tolerance: A delta to add to the bounding box
            coord: (x,y) tuple

         *Returns*:
            Tuple of the bounding box (xmin, ymin, xmax, ymax)

        """

        xmin = coord[0] - tolerance
        ymin = coord[1] - tolerance
        xmax = coord[0] + tolerance
        ymax = coord[1] + tolerance

        return (xmin, ymin, xmax, ymax)

    @staticmethod
    def calculate_adjusted_area(area, cmp_index):
        """Calculate the adjusted area from the area and compactness index

        Args:
            area (float): Area of the polygon
            cmp_index (float): Compactness index of the areea

        Return:
            flot: Adjusted area of the polygon

            """
        return area * (0.75 / cmp_index)

    @staticmethod
    def read_in_file(in_file, geo_content, layer_in=None):
        """
        Read and load the vectors in the input file

        Args:
            in_file (str): Name of the input file (geopackage)
            geo_content (dict): Dictionary containing information to create the spatial database
            layer_in (str): Layer name to read

        Return:
            None

        """

        if layer_in is None:
            # Extract the name of the layers in the file
            geo_content.layer_names = fiona.listlayers(in_file)
        else:
            # Only extract specific layer
            geo_content.layer_names = layer_in

        # extract the spatial feature in the file
        for layer_name in geo_content.layer_names:
            with fiona.open(in_file, 'r', layer=layer_name) as src:
                geo_content.crs = src.crs
                geo_content.driver = src.driver
                geo_content.schemas[layer_name] = src.schema
                geo_content.bounds.append(src.bounds)

                for in_feature in src:
                    geom = in_feature['geometry']
                    if geom['type'] == 'Point':
                        feature = PointSc(geom['coordinates'])
                    elif geom['type'] == 'LineString':
                        feature = LineStringSc(geom['coordinates'])
                    elif geom['type'] == 'Polygon':
                        exterior = geom['coordinates'][0]
                        interiors = geom['coordinates'][1:]
                        feature = PolygonSc(exterior, interiors)
                    else:
                        print("The following geometry type is unsupported: {}".format(geom['type']))
                        feature = None
                    if feature is not None:
                        feature.sb_layer_name = layer_name  # Layer name is the key for the schema
                        feature.sb_properties = in_feature['properties']
                        geo_content.in_features.append(feature)
            src.close()

    @staticmethod
    def write_out_file(out_file, geo_content):
        """
        Write the vectors in the output file

        Args:
            out_file (str): Name of the output file (geopackage)
            geo_content (DataClass): Contains information to create the spatial database

        Return:
            None

        """

        # Loop over each layer and write the content of the file
        for layer_name in geo_content.layer_names:
            with fiona.open(out_file, 'w',
                            driver=geo_content.driver,
                            layer=layer_name,
                            crs=geo_content.crs,
                            schema=geo_content.schemas[layer_name]) as dest:
                out_features = []
                for feature in (feature for feature in geo_content.out_features
                                if feature.sb_layer_name == layer_name):
                    # Transform the Shapely features for fiona writing
                    if feature.geom_type == GenUtil.POINT:
                        coordinates = (feature.x, feature.y)
                        geo_content.out_nbr_points += 1
                    elif feature.geom_type == GenUtil.LINE_STRING:
                        coordinates = list(feature.coords)
                        geo_content.out_nbr_line_strings += 1
                    elif feature.geom_type == GenUtil.POLYGON:
                        exterior = list(feature.exterior.coords)
                        interiors = [list(interior.coords) for interior in feature.interiors]
                        coordinates = [exterior] + interiors
                        geo_content.out_nbr_polygons += 1
                        geo_content.out_nbr_holes += len(interiors)

                    out_feature = {'geometry': {'type': feature.geom_type,
                                                'coordinates': coordinates},
                                   'properties': feature.sb_properties}
                    out_features.append(out_feature)

                dest.writerecords(out_features)

            dest.close()

    @staticmethod
    def write_out_file_append(out_file, geo_content):
        """
        Write the vectors in the output file

        Args:
            out_file (str): Name of the output file (geopackage)
            geo_content (DataClass): Contains information to create the spatial database

        Return:
            None

        """

        line_schema = landmarks_schema = {'geometry': 'LineString',
                                          'properties': OrderedDict([])
                                          }

        # Loop over each layer and write the content of the file
        for layer_name in geo_content.layer_names:
            with fiona.open(out_file, 'w',
                            driver=geo_content.driver,
                            layer=layer_name,
                            crs=geo_content.crs,
                            schema=line_schema) as dest:
                out_features = []
                for feature in (feature for feature in geo_content.out_features
                                if feature.sb_layer_name == layer_name):
                    # Transform the Shapely features for fiona writing
                    if feature.geom_type == GenUtil.POINT:
                        coordinates = (feature.x, feature.y)
                        geo_content.out_nbr_points += 1
                    elif feature.geom_type == GenUtil.LINE_STRING:
                        coordinates = list(feature.coords)
                        geo_content.out_nbr_line_strings += 1
                    elif feature.geom_type == GenUtil.POLYGON:
                        exterior = list(feature.exterior.coords)
                        interiors = [list(interior.coords) for interior in feature.interiors]
                        coordinates = [exterior] + interiors
                        geo_content.out_nbr_polygons += 1
                        geo_content.out_nbr_holes += len(interiors)

                    out_feature = {'geometry': {'type': feature.geom_type,
                                                'coordinates': coordinates},
                                   'properties': feature.sb_properties}
                    out_features.append(out_feature)

                dest.writerecords(out_features)

            dest.close()


class SpatialContainer(object):
    """This class manages the spatial features and a spatial index.

    This class enables the management of spatial features by incorporation
    transparently a spatial index.  The spatial index is an implementation
    of the Rtree open source softawre.  The spatial container offers the following
    main options:
      - add features in the constainer and update the spatial index
      - delete features in the container and update the spatial index
      - update the coordinates of a feature and update the spatial index if needed
      - make spatial queries by bounding box
      - make attributes queries
      - delete the container

    """

    # Class variable that contains the Spatial Container Internal ID
    _sc_id = 0

    def __init__(self):
        """Create an object of type SpatialContainer

        The init will create one container for the feature a dictionary and one
        container for the spatial index (Rtree)

        *Parameters*: None

        *Returns*: *None*

        """

        self._r_tree = Rtree()  # Container for the Rtree
        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def adjust_bbox(self, bounds, delta=GenUtil.ZERO):
        """Adjust the bounding box by increasing by a very small delta

        Parameters:
            bounds: Tuple forming the bounding box (xmin, ymin, wmax, ymax)

        return value:
            altered bounding box (xmin, ymin, wmax, ymax)"""

        xmin, ymin, xmax, ymax = bounds

        xmin -= delta
        ymin -= delta
        xmax += delta
        ymax += delta

        return (xmin, ymin, xmax, ymax)

    def add_feature(self, feature):
        """Adds a feature in the container and update the spatial index with the feature's bound

        To be added in the container a spatial feature must be a MA_Point, MA_LineString or
        MA_Polygon.

        *Parameters*:
            - feature: A spatial feature derives from PointSc, LineStringSc

        *Returns*: *None*

        """

        # Check if the type is valid
        if issubclass(feature.__class__, (PointSc, LineStringSc, PolygonSc)):
            pass
        else:
            raise GeoSimException('Unsupported class: {}'.format(str(feature.__class__)))

        bounds = feature.bounds

        # Adjust the bounding box
        bounds = self.adjust_bbox(bounds)

        # Container unique internal counter
        SpatialContainer._sc_id += 1

        # Add the spatial id to the feature
        feature._sc_id = SpatialContainer._sc_id
        feature._sc_scontainer = self

        # Add the feature in the feature container
        self._features[feature._sc_id] = feature

        # Add the bounding box in the bbox_container
        self._bbox_features[feature._sc_id] = bounds
        self._r_tree.add(feature._sc_id, bounds)

        return

    def add_features(self, features):
        """Adds a list of feature in the spatial container and

        *Parameters*:
            - feature: A spatial feature derived from Point, LineString or Polygon

        *Returns*: *None*

        """

        for feature in features:
            self.add_feature(feature)

        return

    def del_feature(self, feature):
        """Delete the feature in the spatial container and in the RTree.

        If the feature is included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...

        *Parameters*:
            - feature: The feature to delete in the spatial container

        *Returns*:
            None

        """

        ret_value = 0

        # Check if the feature has a container_key
        if hasattr(feature, "_sc_id"):

            if (feature._sc_id in self._features and
                    feature._sc_id in self._bbox_features):

                try:
                    # Retrieve the bounding boxes of this feature
                    bbox = self._bbox_features[feature._sc_id]
                    # Delete the feature from the features and the bbox_features
                    del self._features[feature._sc_id]
                    del self._bbox_features[feature._sc_id]
                    # Delete the different bounds in RTree
                    self._r_tree.delete(feature._sc_id, bbox)
                    # Reset the property _sc_id and _sc_scontainer
                    feature._sc_id = None
                    feature._sc_scontainer = None
                except:
                    raise InternalError("Internal corruption, problem with the container and/or the RTree")
            else:
                raise InternalError("Internal corruption, key {} has disappear...".format(feature._sc_id))

        return ret_value

    def del_features(self, features):
        """Delete a list of features in the spatial container

        If the features are included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...

        *Parameters*:
            - features: list of features to delete

        *Returns*:
            - List of value for each feature to delete.
            - 0 if feature is deleted from the patial container
            - 1 if feature was not included in the spatial container

        Exception InternalError: If the key is not in one of the structure

        """

        ret_value = []

        for feature in features:
            ret_value.append(self.del_feature(feature))

        return ret_value

    def update_spatial_index(self, feature):
        """Update the bounds of the feature in the spatial index

        It will only modify the Rtree spatial index if the bounding
        box of the feature is changed in comparison with the old one.

        *Parameters*:
            - feature: Feature containing the bounds to update

        *Returns*: *None*

        """

        old_bbox = self._bbox_features[feature._sc_id]
        new_bbox = feature.bounds
        old_x_min, old_y_min, old_x_max, old_y_max = old_bbox[0], old_bbox[1], old_bbox[2], old_bbox[3]
        new_x_min, new_y_min, new_x_max, new_y_max = new_bbox[0], new_bbox[1], new_bbox[2], new_bbox[3]

        if old_x_min <= new_x_min and \
                old_y_min <= new_y_min and \
                old_x_max >= new_x_max and \
                old_y_max >= new_y_max:
            # Nothing to do new bounding box is completely included into the old one
            pass
        else:
            # The bounding box has changed
            # Adjust the bounding box
            new_bbox = self.adjust_bbox(new_bbox)
            # Delete The old bounding box in Rtree
            self._r_tree.delete(feature._sc_id, old_bbox)
            # Add the new bounding boxes in Rtree
            self._r_tree.add(feature._sc_id, new_bbox)

            # Save the bounding boxes
            self._bbox_features[feature._sc_id] = new_bbox

        return

    def get_features(self, bounds=None, remove_features=[]):
        """Extract the features from the spatial container.

        According to the parameters the extraction can manage the extraction based on a bounding box using
        the spatial index RTree, some filters to manage extraction based on properties and the possibility
        to remove specific features based on a list of keys

        *Parameters*:
            - bounds: Bounding for the spatial extraction. *None* means all the features
            - remove_keys: List of keys to be removed from the selection

        *Returns*:
            - List of features extracted from spatial container

        """

        tmp_remove_features = []
        for feature in remove_features:
            if isinstance(feature, int):
                tmp_remove_features.append(feature)
            else:
                tmp_remove_features.append(feature._sc_id)

        remove_features = tmp_remove_features

        # Extract the features by bounds if requested
        if bounds is not None:
            # Extract features by bounds
            keys = (key for key in self._r_tree.intersection(bounds) if key not in remove_features)
            features = [self._features[key] for key in keys if key in self._features]
        else:
            features = [feature for feature in self._features.values() if feature not in remove_features]

        # The following code allows to delete feature while iterating over a get_features request
        for feature in features:
            if feature._sc_id is not None:
                yield feature
            else:
                # If the feature has been deleted do not return it
                pass

        return


class SpatialContainerSTRtree(object):
    """This class manages the spatial features and a spatial index for the STRtree.

    The STRtree is using the STRtree of shapely and has the following limitations
    compared to RTree:
       - immutable rtree
       - only able load features one time (all at the same time)
       - no edit after
       - edition is needed an exception is thrown

    """

    # Class variable that contains the Spatial Container Internal ID
    _sc_id = 0

    def __init__(self):
        """Create an object of type SpatialContainer

        The init will create one container for the feature a dictionary and one
        container for the spatial index (Rtree)

        *Parameters*: None

        *Returns*: *None*

        """

        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def adjust_bbox(self, bounds, delta=GenUtil.ZERO):
        """Adjust the bounding box by increasing by a very small delta

        Parameters:
            bounds: Tuple forming the bounding box (xmin, ymin, wmax, ymax)

        return value:
            altered bounding box (xmin, ymin, wmax, ymax)"""

        xmin, ymin, xmax, ymax = bounds

        xmin -= delta
        ymin -= delta
        xmax += delta
        ymax += delta

        return (xmin, ymin, xmax, ymax)

    def add_features(self, features):
        """Adds all the features in the container and update the spatial index with the feature's bound


        *Parameters*:
            - feature: A spatial feature derives from PointSc, LineStringSc

        *Returns*: *None*

        """

        tmp_features = []
        for feature in features:
            # Check if the type is valid
            if issubclass(feature.__class__, (PointSc, LineStringSc, PolygonSc)):
                pass
            else:
                raise GeoSimException('Unsupported class: {}'.format(str(feature.__class__)))

            bounds = feature.bounds

            # Adjust the bounding box
            bounds = self.adjust_bbox(bounds)

            # Container unique internal counter
            SpatialContainer._sc_id += 1

            # Add the spatial id to the feature
            feature._sc_id = SpatialContainer._sc_id
            feature._sc_scontainer = self

            # Add the feature in the feature container
            self._features[feature._sc_id] = feature

            # Add the bounding box in the bbox_container
            self._bbox_features[feature._sc_id] = bounds

            # Transform the feature as its corresponding bounding box... to simulate the Rtree class
            xmin, ymin, xmax, ymax = bounds
            tmp_feature = LineString(((xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)))
            tmp_feature._sc_id = feature._sc_id
            tmp_features.append(tmp_feature)

        # Load all the features at the same time in the shapely rtree
        self._r_tree = STRtree(tmp_features)

        return

    def add_feature(self, features):
        """Adds a list of feature in the spatial container and

        *Parameters*:
            - feature: A spatial feature derived from Point, LineString or Polygon

        *Returns*: *None*

        """

        raise GeoSimException("Cannot add feature with shapely STRtree")

        return

    def del_feature(self, feature):
        """Delete the feature in the spatial container and in the RTree.

        Raise exception cannot delete feature in shapely RTree

        *Parameters*:
            - feature: The feature to delete in the spatial container

        *Returns*:
            None

        """

        raise GeoSimException("Cannot delete feature with shapely STRtree")

    def del_features(self, features):
        """Delete a list of features in the spatial container

        Raise exception cannot delete feature in shapely RTree

        *Parameters*:
            - features: list of features to delete

        *Returns*:
            Exception

        Exception InternalError: If the key is not in one of the structure

        """

        raise GeoSimException("Cannot delete features with shapely STRtree")

        return ret_value

    def update_spatial_index(self, feature):
        """Update the bounds of the feature in the spatial index

        It will only modify the Rtree spatial index if the bounding
        box of the feature is changed in comparison with the old one.

        *Parameters*:
            - feature: Feature containing the bounds to update

        *Returns*: *None*

        """

        old_bbox = self._bbox_features[feature._sc_id]
        new_bbox = feature.bounds
        old_x_min, old_y_min, old_x_max, old_y_max = old_bbox[0], old_bbox[1], old_bbox[2], old_bbox[3]
        new_x_min, new_y_min, new_x_max, new_y_max = new_bbox[0], new_bbox[1], new_bbox[2], new_bbox[3]

        if old_x_min <= new_x_min and \
                old_y_min <= new_y_min and \
                old_x_max >= new_x_max and \
                old_y_max >= new_y_max:
            # Nothing to do new bounding box is completely included into the old one
            pass
        else:
            # The bounding box has changed...
            raise GeoSimException("Cannot change the bounding box of a feature with shapely STRtree")

        return

    def get_features(self, bounds=None, remove_features=[]):
        """Extract the features from the spatial container.

        According to the parameters the extraction can manage the extraction based on a bounding box using
        the spatial index RTree, some filters to manage extraction based on properties and the possibility
        to remove specific features based on a list of keys

        *Parameters*:
            - bounds: Bounding for the spatial extraction. *None* means all the features
            - remove_keys: List of keys to be removed from the selection

        *Returns*:
            - List of features extracted from spatial container

        """

        tmp_remove_features = []
        for feature in remove_features:
            if isinstance(feature, int):
                tmp_remove_features.append(feature)
            else:
                tmp_remove_features.append(feature._sc_id)

        remove_features = tmp_remove_features

        # Extract the features by bounds if requested
        if (bounds != None):
            # Extract features by bounds
            keys = (key for key in self._r_tree.intersection(bounds) if key not in remove_features)
            features = [self._features[key] for key in keys if key in self._features]
        else:
            features = [feature for feature in self._features.values() if feature not in remove_features]

        # The following code allows to delete feature while iterating over a get_features request
        for feature in features:
            if feature._sc_id is not None:
                yield feature
            else:
                # If the feature has been deleted do not return it
                pass

        return


class ChordalAxis(object):

    def __init__(self, lst_triangle, search_tolerance=GenUtil.ZERO):
        """Constuctor of the ChordalAxis class

        *Parameters*:
            - lst_triangle: List of LineString triangle
            - search_tolerance: float for the zero value approximation

        *Returns*:
            -

                """

        self.search_tolerance = search_tolerance

        self._validate_triangles(lst_triangle)

        # Transform the Triangle LineString into LineStringSc to be loaded in SpatialContainer
        for i, triangle in enumerate(lst_triangle):
            triangle = LineStringSc(triangle.coords)
            triangle._id = i
            lst_triangle[i] = triangle

        self._build_topology(lst_triangle)

        self.triangle_clusters = self._build_clusters()

        return

    def _validate_triangles(self, lst_triangle):
        """Validate triangles geometry

        *Parameters*:
            - lst_triangle: List of LineString of triangle

        *Returns*:
            - None

        """

        if len(lst_triangle) >= 1:
            triangle_type = lst_triangle[0].geom_type
        triangle_valid = True
        for i, triangle in enumerate(lst_triangle):

            # Triangle are LineString
            if triangle.geom_type == GenUtil.LINE_STRING:
                coords = list(triangle.coords)
                # Check number of vertice
                if len(coords) != 4:
                    print("Triangle does not contain exactly 4 vertices: {0}".format(coords[0]))
                    triangle_valid = False
                # Check if all geometry are identical
                if triangle.geom_type != triangle_type:
                    print("Triangle has mixed geometry type: {0}".format(coords[0]))
                    triangle_valid = False
                # Check if the line is closed
                if Point(coords[0]).distance(Point(coords[3])) >= self.search_tolerance:
                    print("Triangle is not closed: {0}".format(coords[0]))
                    triangle_valid = False

            # Triangle are polygons
            if triangle.geom_type == GenUtil.POLYGON:
                coords = list(triangle.exterior.coords)
                # Validate trianglehas 4 vertice
                if len(coords) != 4:
                    print("Triangle does not contain exactly 4 vertices: {0}".format(coords[0]))
                    triangle_valid = False
                    # Check if all geometry are identical
                if triangle.geom_type != triangle_type:
                    print("Triangle has mixed geometry type: {0}".format(coords[0]))
                    triangle_valid = False

                # Transform the polygon into a LineString
                if triangle_valid:
                    triangles[i] = LineString(coords)

            if not triangle_valid:
                # There are one or more errors in the triangles
                raise GeoSimException("Error in the triangles... cannot process them...")

            return

    def _find_adjacent_triangle(self, in_hand_triangle, mid_pnt_side):
        """Test if the parameter is iterable; if not make it iterable by creating a tuple of one element

        *Parameters*:
            - in_hand_triangle: Triangle (LineStringSc) to check for adjacency
            - mid_pnt_side: Point located on one of the edge of the triangle used as search point for adjacency

        *Returns*:
            - LineStringSc of the adjacent triangle or None (if no triangle)

        """

        # Find adjacent triangles
        triangles = self.s_container.get_features(bounds=mid_pnt_side.bounds,remove_features=[in_hand_triangle])

        min_distance = self.search_tolerance
        target_triangle = None
        nbr_near_zero = 0
        for triangle in triangles:
            if triangle._id != in_hand_triangle._id:
                distance = mid_pnt_side.distance(triangle)
                if distance < self.search_tolerance:
                    nbr_near_zero += 1
                    if distance < min_distance:
                        min_distance = distance
                        target_triangle = triangle
            else:
                # Do not process same triangle
                pass

        if nbr_near_zero >= 2:
            xy = (mid_pnt_side.x, mid_pnt_side.y)
            print("***Warning*** Triangles are to small: {0} Try to simplify them (e.g. Douglas Peucker)".format(xy))

        return target_triangle

    def _set_attributes(self, triangle):
        """Sets some attribute of the triangle

        *Parameters*:
            - triangle: LineStringSC to calculate attributes

        *Returns*:
            - None

        """

        # Add attributes to the triangle
        triangle._centre_lines = []  # List to store center line
        triangle._adjacent_sides = []  # counter of the number of adjacent side (max=3)

        # Calculate the mid point of each side of the triangle
        coords = list(triangle.coords)
        mid_pnt_side_0 = LineString([coords[0], coords[1]]).interpolate(0.5, normalized=True)
        mid_pnt_side_1 = LineString((coords[1], coords[2])).interpolate(0.5, normalized=True)
        mid_pnt_side_2 = LineString((coords[2], coords[0])).interpolate(0.5, normalized=True)
        triangle._mid_pnt_side = [mid_pnt_side_0, mid_pnt_side_1, mid_pnt_side_2]

        # Finsd adjacent triangle on each side
        adjacent_side_0 = self._find_adjacent_triangle(triangle, mid_pnt_side_0)
        adjacent_side_1 = self._find_adjacent_triangle(triangle, mid_pnt_side_1)
        adjacent_side_2 = self._find_adjacent_triangle(triangle, mid_pnt_side_2)
        triangle._adjacent_sides_ref = [adjacent_side_0, adjacent_side_1, adjacent_side_2]

        # Set a list of 3 items: 0: No triangle adjacent; 1: Yes triangle adjacent
        for adjacent_sides_ref in triangle._adjacent_sides_ref:
            if adjacent_sides_ref is not None:
                triangle._adjacent_sides.append(1)  # There is a triangle on this side
            else:
                triangle._adjacent_sides.append(0)  # There is no triangle on this side
        triangle._type = sum(triangle._adjacent_sides)  # Number of side adjacent to another triangle (max: 3)

        return

    def _build_topology(self, lst_triangle):
        """Build the topology

        Determine for each triangle if there is an adjacent triangle

        *Parameters*:
            - lst_triangles: list of LineStringSc

        *Returns*:
            - None

        """

        # Create spatial container
        self.s_container = SpatialContainer()

        # Load triangles
        self.s_container.add_features(lst_triangle)

        # Loop to set ssome attributes on each triangle
        for triangle in self.s_container.get_features():
            self._set_attributes(triangle)

        return

    def _build_clusters(self):
        """Build the clusters of triangle

        One cluster of triangles are all the triangles that have a common edge. One cluster of polygon is equivalent
        to the area of one polygon (including the holes). The set of clusters represent all the polygons

        *Parameters*:
            - None

        *Returns*:
            - Clusters: list of list of triangles

        """

        dict_triangles = {}
        clusters = []
        for triangle in self.s_container.get_features():
             dict_triangles[triangle._id] = triangle

        # Loop unitl all the triangles are processed
        while len(dict_triangles) >= 1:
            seed_triangle = next(iter(dict_triangles.values()))
            cluster = self._build_one_cluster(dict_triangles, seed_triangle)
            clusters.append(cluster)

        return clusters

    def _build_one_cluster(self, dict_triangles, seed_triangle):
        """Identify all the triangle that shared an edge (triangle part of one polygon)

        This method is simulating recursivity using a stack

        *Parameters*:
            - dict_triangle: Dictionary of the triangle
            - seed_triangle: Random starting point to form a cluster

        *Returns*:
            - List of triangle forming a cluster

        """

        # Create cluster list to accumulate triangle in cluster
        cluster = []

        # Create the stack to simulate recursivity
        stack = []

        # Initialize the stack with the seed triangle
        stack.append(seed_triangle)

        # Loop over the stack until no more triangle to process
        while stack:
            # Fetch the next triangle
            triangle = stack.pop()
            # Add the tiangle in the cluster
            cluster.append(triangle)

            if triangle._id in dict_triangles:
                # Remove from the triangle from the dictionary
                del dict_triangles[triangle._id]
                # Process the adjacent sides
                for adjacent_side_ref in (triangle._adjacent_sides_ref):
                    if adjacent_side_ref is not None:
                        if adjacent_side_ref._id in dict_triangles:
                            stack.append(adjacent_side_ref) # Simulate recursivity
                        else:
                            # triangle already processed
                            pass
                    else:
                        # No triangle to process
                        pass
            else:
                # Triangle alrerady processed
                pass

        return cluster

    def _create_centre_line(self, triangle):
        """Calculates and extract the center line of one triangle

        The center line depends of the type of triangle
            Terminal triangle: From the side of the adjacent edge to the oposite angle
            Sleeve triangle: Joining the mid side of the internal side
            Junction triangle: Calculate the centroid and create 3 lines from the centroid to the mid of each side

        *parameters*:
            - triangle: LineStringSc triangle

        *Returns*:
            - None

        """

        coords = list(triangle.coords)

        # Process each case depending on the number of internal side of the triangle
        if triangle._type == 0:
            # Degenerated polygon with one triangle no skeleton line added
            pass

        if triangle._type == 1:
            # Terminal triangle add line from the extremity of the triangle up to mid opposite side
            if triangle._adjacent_sides[0] == 1:
                coords_line = [coords[2], triangle._mid_pnt_side[0]]
            if triangle._adjacent_sides[1] == 1:
                coords_line = [coords[0], triangle._mid_pnt_side[1]]
            if triangle._adjacent_sides[2] == 1:
                coords_line = [coords[1], triangle._mid_pnt_side[2]]

            triangle._centre_lines.append(LineString(coords_line))

        if triangle._type == 2:
            # Sleeve triangle skeleton added between the mid point of side adjacent to another triangle
            mid_pnt = []
            for i,side in enumerate(triangle._adjacent_sides):
                if side == 1:
                    mid_pnt.append(triangle._mid_pnt_side[i])
            triangle._centre_lines.append(LineString([mid_pnt[0], mid_pnt[1]]))

        if triangle._type == 3:
            # Junction triangle T type skeleton added.
            centroid_x = (coords[0][0] + coords[1][0] + coords[2][0]) / 3.
            centroid_y = (coords[0][1] + coords[1][1] + coords[2][1]) / 3.
            centroid = [centroid_x, centroid_y]

            for mid_side_pnt in triangle._mid_pnt_side:
                triangle._centre_lines.append(LineString([centroid, mid_side_pnt]))

        return

    def get_skeleton(self):
        """extract the ceneter line of each triangle merged them and create a list of LineString feature

                This method is simulating recursivity using a stack

        *Parameters*:
            - None

        *Returns*:
            - None

        """

        merged_centre_lines = []

        # Process each cluster (polygon) independently
        for triangle_cluster in self.triangle_clusters:
            centre_lines = []
            # Process each triangle of one cluster
            for triangle in triangle_cluster:
                self._create_centre_line(triangle)
                centre_lines += triangle._centre_lines

            merge_centre_line = linemerge(centre_lines)
            merged_centre_line = GenUtil.make_iterable(merge_centre_line)
            merged_centre_lines += merged_centre_line

        return merged_centre_lines
    

class GeoSimException(Exception):
    """
    This is the base exception class for genmetal algorithms
    """

    def __init__(self, *arguments, **keywords):
        Exception.__init__(self, *arguments, **keywords)


class InternalError(GeoSimException):
    """
    This exception is raised when an internal error as occurred
    """

    def __init__(self, *param_names):
        """
        Initialise an Invalid Geometry Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GeoSimException.__init__(self, *param_names)
