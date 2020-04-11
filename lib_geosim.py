#!/usr/bin/env python
# -=- encoding: utf-8 -=-

"""
General classes and utilities needed for the GeoSim.

"""

import sys
from math import atan, degrees
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import linemerge
from shapely.ops import unary_union
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

    # Define the type of Triangle
    ISOLATED = 0  # No side touch another triangl
    TERMINAL = 1  # One side touche another triangle
    SLEEVE = 2    # Two sides touch a triangle
    JUNCTION = 3  # All sides touch a triangle

    # Define the type of action(type_action) for creating the centre line
    NONE = 0  # No special action is needed when creating the center line
    T_EDIT_CENTRE_LINE = 1  # A correction is needed for a T Junction (only used when correcting the centre line)
    X_EDIT_CENTRE_LINE = 2  # A correction is needed for a X Junction (only used when correcting the centre line)
    NO_CENTRE_LINE = 3  # No centre line is needed for this triangle (only used when correcting the centre line)

    ANGLE_T_JUNCTION = 45.  # Delta used to test if 2 branches or contiguous
    SEARCH_TOLERANCE = None

    def __init__(self, lst_triangle, search_tolerance=GenUtil.ZERO):
        """Constuctor of the ChordalAxis class

        *Parameters*:
            - lst_triangle: List of LineString triangle
            - search_tolerance: float for the zero value approximation

        *Returns*:
            -

        """

        ChordalAxis.SEARCH_TOLERANCE = search_tolerance

        self._validate_triangles(lst_triangle)

        # Transform the Triangle LineString into _TriangleSc to be loaded in SpatialContainer
        for i, triangle in enumerate(lst_triangle):
            triangle = _TriangleSc(triangle.coords)
            lst_triangle[i] = triangle

        # Create spatial container
        self.s_container = SpatialContainer()

        # Load triangles
        self.s_container.add_features(lst_triangle)

        # Load some class viriables
        _TriangleSc.s_container = self.s_container

        # Build the cluster (group of triangles part of a polygon)
        self.triangle_clusters = self._build_clusters()

        self.nbr_polygons = len(self.triangle_clusters)
        self.nbr_triangles = len(lst_triangle)
        self.nbr_lines_pruned = 0
        self.nbr_iteration = 0
        self.nbr_t_junction = 0
        self.nbr_x_junction = 0

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
                if Point(coords[0]).distance(Point(coords[3])) >= ChordalAxis.SEARCH_TOLERANCE:
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
                    lst_triangle[i] = LineString(coords)

            if not triangle_valid:
                # There are one or more errors in the triangles
                raise GeoSimException("Error in the triangles... cannot process them...")

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
             dict_triangles[triangle.id] = triangle

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

            if triangle.id in dict_triangles:
                # Remove from the triangle from the dictionary
                del dict_triangles[triangle.id]
                # Process the adjacent sides
                for adjacent_side_ref in (triangle.adjacent_sides_ref):
                    if adjacent_side_ref is not None:
                        if adjacent_side_ref.id in dict_triangles:
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
                centre_lines += triangle.centre_line

            merge_centre_line = linemerge(centre_lines)
            merged_centre_line = GenUtil.make_iterable(merge_centre_line)
            merged_centre_lines += merged_centre_line

        return merged_centre_lines

    def correct_skeleton(self):

        # Prune the small branch from the Junction triangle until there are no more small branch to prune
        while True:
            self.nbr_iteration += 1
            nbr_pruned = 0
            for triangle in self.s_container.get_features():
                if triangle.type == ChordalAxis.JUNCTION:
                    nbr_pruned += self.prune_junction(triangle)
            if nbr_pruned == 0:
                break
            else:
                self.nbr_lines_pruned += nbr_pruned

        # Correct the X junction to join adjacent junction and remove the line joining the 2 junctions
        for triangle in self.s_container.get_features():
            if triangle.type == ChordalAxis.JUNCTION and triangle.sub_type == ChordalAxis.NONE:
                sides_x_junction = self.adjust_x_junction(triangle)
                if sides_x_junction is not None:
                    adjacent_triangles = sides_x_junction[0]
                    mid_sides_pnt = sides_x_junction[1]
                    centroid = sides_x_junction[2]

                    self.nbr_x_junction += 1
                    triangle.sub_type = ChordalAxis.X_EDIT_CENTRE_LINE
                    triangle.class_x_mid_sides_pnt = mid_sides_pnt
                    triangle.class_x_centroid = centroid
                    triangle.centre_line = None
                    for adjacent_triangle in adjacent_triangles:
                        adjacent_triangle.sub_type = ChordalAxis.NO_CENTRE_LINE
                        adjacent_triangle.centre_line = None
                else:
                    # Not a valid X junction. Nothing to do
                    pass

        # Correct the T junction to form a straight line if 2 brancj have the same orientation
        for triangle in self.s_container.get_features():
            if triangle.type == ChordalAxis.JUNCTION and triangle.sub_type == ChordalAxis.NONE:
                sides_t_junction = self.adjust_t_junction(triangle)
                if sides_t_junction is not None:
                    self.nbr_t_junction += 1
                    triangle.sub_type = ChordalAxis.T_EDIT_CENTRE_LINE
                    triangle.junction_side_a = sides_t_junction[0]
                    triangle.junction_side_b = sides_t_junction[1]
                    triangle.centre_line = None


    def adjust_t_junction(self, junction_triangle):
        """This function corrects T junction"""

        branches = []
        sides_t_junction = None
        for next_triangle in junction_triangle.adjacent_sides_ref:
            if next_triangle.type == ChordalAxis.SLEEVE:
                branch = Branch(junction_triangle, next_triangle)
                branches.append(branch)

        branch_angle = [branch.angle for branch in branches] # Extract the angles of the branch
        if len(branches) == 3 and None not in [branch_angle]:
            sides_t_junction = None
            angle_max = ChordalAxis.ANGLE_T_JUNCTION
            for i,j in [(0,1),(1,2),(2,0)]:
                delta_angle = abs(180. - abs(branch_angle[i] - branch_angle[j]))
                if delta_angle < angle_max:
                    angle_max = delta_angle
                    sides_t_junction = [i, j]
            else:
                # No correction done
                pass
        else:
            # No T junction to process
            pass

        return sides_t_junction

    def adjust_x_junction(self, current_junction):

        side_x_junction = None
        for adjacent_junction in current_junction.adjacent_sides_ref:
            branch = Branch(current_junction, adjacent_junction)
            last_triangle = branch.triangle_in_branch[-1]  # Extract the last triangle of the branch

            if last_triangle.type == ChordalAxis.JUNCTION and \
               last_triangle.sub_type == ChordalAxis.NONE and \
               branch.length < current_junction.width:

                # Merge the triangle to form only one polygon
                line_triangles = [current_junction] + branch.triangle_in_branch
                pol_triangles = [Polygon(line.coords) for line in line_triangles]
                merged_pol = unary_union(pol_triangles)
                if merged_pol.geom_type == GenUtil.POLYGON:  # Merged the triangle to form only one polygon
                    centroid = merged_pol.centroid
                    merged_line = LineString(merged_pol.exterior.coords)

                    # Detect which mid side point we must keep
                    mid_pnt_sides = current_junction.mid_pnt_sides + last_triangle.mid_pnt_sides
                    new_mid_pnt_sides = []
                    for mid_pnt_side in mid_pnt_sides:
                        if mid_pnt_side.distance(merged_line) < ChordalAxis.SEARCH_TOLERANCE:
                            new_mid_pnt_sides.append(mid_pnt_side)

                    if self.validate_junction_x(merged_pol, centroid, new_mid_pnt_sides):
                        side_x_junction = [branch.triangle_in_branch, new_mid_pnt_sides, centroid]
                    break
                else:
                    # Invalid merging nothing to do
                    pass
            else:
                # Not a triangle to process for X junction
                pass

        return side_x_junction

    def validate_junction_x(self, merged_pol, centroid, new_mid_pnt_sides):

        buf_merged_pol = merged_pol.buffer(.01)
        print ("Validate x junction add widht function")
        status = True
        for mid_pnt_side in new_mid_pnt_sides:
            line = LineString((centroid, mid_pnt_side))
            if line.crosses(buf_merged_pol):
                valid= False
                break

        return status


    def prune_junction(self, junction_triangle):
        """This function prune a junction triangle of the branches that are below a certain tolerance"""

        branches = []

        for next_triangle in junction_triangle.adjacent_sides_ref:
            branch = Branch(junction_triangle, next_triangle)
            # Only keep small TERMINAL branches
            if branch.last_triangle_type == ChordalAxis.TERMINAL and branch.length <= junction_triangle.width:
                branches.append(branch)

        if len(branches) == 3:
            # The three branches of the junction Triangle are below the tolerance
            # Remove the  branch with the smallest tolerance
            max_length = sys.float_info.max
            for branch in branches:
                if branch.length < max_length:
                    del_branches = [branch]
                    max_length = branch.length

        elif len(branches) == 2:
            # Two branches are below the tolerance
            if branches[0].length < branches[1].length:
                branch_0 = branches[0]
                branch_1 = branches[1]
            else:
                branch_0 = branches[1]
                branch_1 = branches[0]
            if branch_0.length < .3 *branch_1.length:
                del_branches = [branch_0]
            else:
                del_branches = [branch_0, branch_1]
        elif len(branches) == 1:
            del_branches = [branches[0]]
        else:
            del_branches = []

        if len(del_branches) >= 1:
            for branch in del_branches:
                # Loop over each branch
                for triangle in branch.triangle_in_branch:
                    if triangle.type in (ChordalAxis.SLEEVE, ChordalAxis.TERMINAL):
                        # Make all the triangle forming the branch has terminal (isolated) terminal
                        triangle.set_as_isolated()

        return len(del_branches)


class _TriangleSc(LineStringSc):

    """LineString specialization to be included in the SpatialContainer"""

    id = 0

    def __init__(self, coords):
        super().__init__(coords)

        # Add unique id to each Triangle
        self.id = _TriangleSc.id

        # Attribute for Junction specialization
        self.sub_type = ChordalAxis.NONE
        _TriangleSc.id += 1


    @property
    def mid_pnt_sides(self):

        try:
            return self._mid_pnt_sides
        except AttributeError:
            # Calculate the mid point of each side of the triangle
            coords = list(self.coords)
            mid_pnt_side_0 = LineString([coords[0], coords[1]]).interpolate(0.5, normalized=True)
            mid_pnt_side_1 = LineString((coords[1], coords[2])).interpolate(0.5, normalized=True)
            mid_pnt_side_2 = LineString((coords[2], coords[0])).interpolate(0.5, normalized=True)
            self._mid_pnt_sides = [mid_pnt_side_0, mid_pnt_side_1, mid_pnt_side_2]
            return self._mid_pnt_sides

    @property
    def type(self):
        try:
            return self._type
        except AttributeError:
            self._type = 0
            for adjacent_side_ref in self.adjacent_sides_ref:
                if adjacent_side_ref != None:
                    self._type += 1

            return self._type

    @property
    def width(self):
        try:
            return self._width
        except AttributeError:
            lst_length = [line.length for line in self.centre_line]
            max_length = max(lst_length)
            self._width = max_length*2.

            return self._width

    @property
    def adjacent_sides_ref(self):
        try:
            return self._adjacent_sides_ref
        except AttributeError:
            self._adjacent_sides_ref = self.get_adjacent_sides_ref()
            return self._adjacent_sides_ref

    @adjacent_sides_ref.setter
    def adjacent_sides_ref(self, value):
        self._adjacent_sides_ref = value

        # Delete different attribute so they are recalculated
        try:
            del self._centre_line
        except AttributeError:
            pass
        try:
           del self._width
        except AttributeError:
            pass
        try:
            del self._type
        except AttributeError:
            pass

    @property
    def centre_line(self):
        try:
            return self._centre_line
        except AttributeError:
            self._centre_line = self._create_centre_line()

            return self._centre_line

    @centre_line.setter
    def centre_line(self, value):
        if value is None:
            try:
                del self._centre_line
            except AttributeError:
                pass
        else:
            self._centre_line = value

    def _create_centre_line(self):
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

        centre_line = []
        coords = list(self.coords)

        # Process each case depending on the number of internal side of the triangle
        if self.type == ChordalAxis.ISOLATED:
            # Degenerated polygon with one triangle no skeleton line added
            pass

        if self._type == ChordalAxis.TERMINAL:
            # Terminal triangle add line from the extremity of the triangle up to mid opposite side
            if self.adjacent_sides_ref[0] != None:
                coords_line = [coords[2], self.mid_pnt_sides[0]]
            if self.adjacent_sides_ref[1] != None:
                coords_line = [coords[0], self.mid_pnt_sides[1]]
            if self.adjacent_sides_ref[2] != None:
                coords_line = [coords[1], self.mid_pnt_sides[2]]

            centre_line.append(LineString(coords_line))

        if self.type == ChordalAxis.SLEEVE:
            if self.sub_type == ChordalAxis.NONE:
                # Sleeve triangle skeleton added between the mid point of side adjacent to another triangle
                mid_pnt = []
                for i,adjacent_side_ref in enumerate(self._adjacent_sides_ref):
                    if adjacent_side_ref != None:
                        mid_pnt.append(self.mid_pnt_sides[i])
                centre_line.append(LineString([mid_pnt[0], mid_pnt[1]]))
            elif self.sub_type == ChordalAxis.NO_CENTRE_LINE:
                # No center line to create
                pass

        if self.type == ChordalAxis.JUNCTION:
            if self.sub_type == ChordalAxis.NONE:
                # Regular triangle T type. Centroid is the centroid of the triangle
                pnt_x = (coords[0][0] + coords[1][0] + coords[2][0]) / 3.
                pnt_y = (coords[0][1] + coords[1][1] + coords[2][1]) / 3.
                centroid = [pnt_x, pnt_y]
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    centre_line.append(LineString([centroid, mid_side_pnt]))

            elif self.sub_type == ChordalAxis.T_EDIT_CENTRE_LINE:
                # Corrected triangle T. Centroid is the middle point between the 2 aligned branches
                pnt0 = self.mid_pnt_sides[self.junction_side_a]
                pnt1 = self.mid_pnt_sides[self.junction_side_b]
                pnt = LineString([(pnt0.x,pnt0.y), (pnt1.x,pnt1.y)]).interpolate(0.5, normalized=True)
                centroid = [pnt.x, pnt.y]
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    centre_line.append(LineString([centroid, mid_side_pnt]))

            elif self.sub_type == ChordalAxis.X_EDIT_CENTRE_LINE:
                centroid = (self.class_x_centroid.x, self.class_x_centroid.y)
                #  create the center line
                for mid_pnt_side in self.class_x_mid_sides_pnt:
                    centre_line.append(LineString ([centroid, mid_pnt_side]) )

            elif self.sub_type == ChordalAxis.NO_CENTRE_LINE:
                # No line to create.  The lines are created by the junction CLASS_X_PRIMARY
                pass

        return centre_line

    def set_as_isolated(self):
        """Set the attribute of the triangle to simulate a terminal triangle"""

        # To simulate an Isolated triangle the current triangle must not reference any adjacent triangle and
        # all adjacent triangle must not reference the current triangle
        # This is the case of a double reference so both reference must be put at None
        for triangle_ref in self.adjacent_sides_ref:
            if triangle_ref != None:
                new_adjacent_sides_ref = triangle_ref.adjacent_sides_ref
                for j in range(3):
                    if new_adjacent_sides_ref[j] is not None and \
                       new_adjacent_sides_ref[j].id == self.id:
                        # Set the first referecne to None
                        new_adjacent_sides_ref[j] = None
                triangle_ref.adjacent_sides_ref = new_adjacent_sides_ref

        # Set the second reference to None
        self.adjacent_sides_ref = [None, None, None]

    def get_adjacent_sides_ref(self):
        """Test if the parameter is iterable; if not make it iterable by creating a tuple of one element

        *Parameters*:
            - in_hand_triangle: Triangle (LineStringSc) to check for adjacency
            - mid_pnt_side: Point located on one of the edge of the triangle used as search point for adjacency

        *Returns*:
            - LineStringSc of the adjacent triangle or None (if no triangle)

        """

        adjacent_sides_ref = []

        # Loop of each side (each mid_pnt_side) to find adjacent triangle
        for mid_pnt_side in self.mid_pnt_sides:

            # Find all potential adjacent triangles
            triangles = _TriangleSc.s_container.get_features(bounds=mid_pnt_side.bounds,remove_features=[self])

            # Find the closest triangle
            min_distance = sys.float_info.max
            nbr_near_zero = 0
            target_triangle = None
            for triangle in triangles:
                distance = mid_pnt_side.distance(triangle)
                if distance < ChordalAxis.SEARCH_TOLERANCE:
                    nbr_near_zero += 1
                    if distance < min_distance:
                        min_distance = distance
                        target_triangle = triangle

            adjacent_sides_ref.append(target_triangle)

            if nbr_near_zero >= 2:
                xy = (mid_pnt_side.x, mid_pnt_side.y)
                print("***Warning*** Triangles may be to small: {0} Try to simplify them (e.g. Douglas Peucker)".format(xy))

        return adjacent_sides_ref


class Branch:
    """Manage one branch 
    """

    def __init__(self, current_triangle, next_triangle):

        self.current_triangle = current_triangle
        self.triangle_in_branch = []
        self.length = 0.
        max_length = current_triangle.width * 3.

        while True:
            self.triangle_in_branch.append(next_triangle)
            if next_triangle.type in (ChordalAxis.SLEEVE, ChordalAxis.TERMINAL):
                self.length += next_triangle.centre_line[0].length  # Add the length
                if next_triangle.type == ChordalAxis.TERMINAL:
                    # This is the end
                    break
            else:
                # It's a junction Triangle
                break

            # Loop for the next adjacent triangle
            if self.length < max_length:
                adjacents = [adjacent for adjacent in next_triangle.adjacent_sides_ref if adjacent is not None]
                if adjacents[0].id == current_triangle.id:
                    current_triangle, next_triangle = next_triangle, adjacents[1]
                else:
                    current_triangle, next_triangle = next_triangle, adjacents[0]
            else:
                # End of looping reached
                break

        # Extract the type of the last triangle in the branch
        self.last_triangle_type = self.triangle_in_branch[-1].type

        return


        # if next_triangle.type == ChordalAxis.JUNCTION:
        #     self.last_triangle = ChordalAxis.JUNCTION
        #     self.triangle_in_branch.append(next_triangle)
        # elif next_triangle.type == ChordalAxis.TERMINAL:
        #     self.last_triangle = ChordalAxis.TERMINAL
        #     self.length = next_triangle.centre_line[0].length
        #     self.triangle_in_branch.append(next_triangle)
        # else:
        #     while next_triangle.type == ChordalAxis.SLEEVE or next_triangle.type == ChordalAxis.TERMINAL:
        #         self.triangle_in_branch.append(next_triangle)  # Add the next triangle in the list
        #         self.length += next_triangle.centre_line[0].length  # Add the length
        #         self.last_triangle = next_triangle.type
        #         if next_triangle.type == ChordalAxis.SLEEVE and self.length < max_length:
        #             # Get the next triangle
        #             adjacents = [adjacent for adjacent in next_triangle.adjacent_sides_ref if adjacent is not None]
        #             if adjacents[0].id == current_triangle.id:
        #                 current_triangle, next_triangle = next_triangle, adjacents[1]
        #             else:
        #                 current_triangle, next_triangle = next_triangle, adjacents[0]
        #         else:
        #             # End of logping reached
        #             break


    @property
    def angle(self):
        try:
            return self._angle
        except AttributeError:
            # Extract the skeleton line from the branch
            lines = []
            for triangle in self.triangle_in_branch:
                if triangle.type in [ChordalAxis.SLEEVE, ChordalAxis.TERMINAL]:
                    lines += triangle.centre_line

            # Merge the lines to form one line
            merged_line = linemerge(lines)
            if merged_line.geom_type == GenUtil.LINE_STRING:
                # Extract the angle formed by the first and last coordinate
                x0, y0 = merged_line.coords[0][0], merged_line.coords[0][-1]
                x1, y1 = merged_line.coords[-1][0], merged_line.coords[-1][-1]

                # Checked that the first coordinate of the line is located on the triangle
                if self.current_triangle.distance(Point(x0,y0)) < ChordalAxis.SEARCH_TOLERANCE:
                    # The line is well oriented
                    pass
                else:
                    # Invert the orientation of the line
                    x0, y0, x1, y1 = x1, y1, x0, y0

                delta_y = y1 - y0
                delta_x = x1 - x0
                if abs(delta_x) <= ChordalAxis.SEARCH_TOLERANCE:  # Avoid division by zero
                    delta_x = ChordalAxis.SEARCH_TOLERANCE
                self._angle = degrees(atan(delta_y / delta_x))

                # Apply a correction for the quadrant; in order to obtain an angle betweer 0..360
                if delta_x >= 0 and delta_y >= 0:
                    # Quadrant 1 nothing to do
                    pass
                elif delta_x < 0 and delta_y >= 0:
                    # Quadrant 2
                    self._angle += 180.
                elif delta_x < 0 and delta_y < 0:
                    # Quadrant 3
                    self._angle += 180.
                else:
                    # delta_x > 0 anf delta_y < 0
                    # Quandrant 4
                    self._angle += 360.
            else:
                # Unknown error... angle is discarded
                self._angle = None

            return self._angle


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
