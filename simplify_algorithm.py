# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module
# pylint: disable=too-many-lines
# pylint: disable=useless-return
# pylint: disable=too-few-public-methods

# /***************************************************************************
# simplify_algorithm.py
# ----------
# Date                 : April 2021
# copyright            : (C) 2020 by Natural Resources Canada
# email                : daniel.pilon@canada.ca
#
#  ***************************************************************************/
#
# /***************************************************************************
#  *                                                                         *
#  *   This program is free software; you can redistribute it and/or modify  *
#  *   it under the terms of the GNU General Public License as published by  *
#  *   the Free Software Foundation; either version 2 of the License, or     *
#  *   (at your option) any later version.                                   *
#  *                                                                         *
#  ***************************************************************************/

"""
QGIS Plugin for Bend reduction
"""


import os
import inspect
import math
from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.core import (QgsProcessing, QgsProcessingAlgorithm, QgsProcessingParameterDistance,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean, QgsFeatureSink, QgsFeatureRequest, QgsPoint,
                       QgsPointXY, QgsLineString, QgsPolygon, QgsWkbTypes, QgsGeometry,
                       QgsGeometryUtils, QgsRectangle, QgsProcessingException, QgsMultiPolygon)
import processing
from .geo_sim_util import Epsilon, GsCollection, GsFeature, GsPolygon, GsLineString, GsPoint


class SimplifyAlgorithm(QgsProcessingAlgorithm):
    """Main class defining the Reduce Bend as a QGIS processing algorithm.
    """

    def tr(self, string):  # pylint: disable=no-self-use
        """Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):  # pylint: disable=no-self-use
        """Returns a new copy of the algorithm.
        """
        return SimplifyAlgorithm()

    def name(self):  # pylint: disable=no-self-use
        """Returns the unique algorithm name.
        """
        return 'simplify'

    def displayName(self):  # pylint: disable=no-self-use
        """Returns the translated algorithm name.
        """
        return self.tr('Simplify')

    def group(self):
        """Returns the name of the group this algorithm belongs to.
        """
        return self.tr(self.groupId())

    def groupId(self):  # pylint: disable=no-self-use
        """Returns the unique ID of the group this algorithm belongs to.
        """
        return ''

    def shortHelpString(self):
        """Returns a localised short help string for the algorithm.
        """
        help_str = """
    Simplify is a geospatial simplification and generalization tool for lines and polygons...

    """

        return self.tr(help_str)

    def icon(self):  # pylint: disable=no-self-use
        """Define the logo of the algorithm.
        """

        cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]
        icon = QIcon(os.path.join(os.path.join(cmd_folder, 'logo.png')))
        return icon

    def initAlgorithm(self, config=None):  # pylint: disable=unused-argument
        """Define the inputs and outputs of the algorithm.
        """

        # 'INPUT' is the recommended name for the main input parameter.
        self.addParameter(QgsProcessingParameterFeatureSource(
                          'INPUT',
                          self.tr('Input layer'),
                          types=[QgsProcessing.TypeVectorAnyGeometry]))

        # 'TOLERANCE' to be used bor bend reduction
        self.addParameter(QgsProcessingParameterDistance(
                          'TOLERANCE',
                          self.tr('Diameter tolerance'),
                          defaultValue=0.0,
                          parentParameterName='INPUT'))  # Make distance units match the INPUT layer units

        # 'VERBOSE' mode for more output information
        self.addParameter(QgsProcessingParameterBoolean(
            'VERBOSE',
            self.tr('Verbose'),
            defaultValue=False))

        # 'OUTPUT' for the results
        self.addParameter(QgsProcessingParameterFeatureSink(
                          'OUTPUT',
                          self.tr('Reduced bend')))

    def processAlgorithm(self, parameters, context, feedback):
        """Main method that extract parameters and call ReduceBend algorithm.
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        # Extract parameter
        source_in = self.parameterAsSource(parameters, "INPUT", context)
        tolerance = self.parameterAsDouble(parameters, "TOLERANCE", context)
        validate_structure = self.parameterAsBool(parameters, "VALIDATE_STRUCTURE", context)
        verbose = self.parameterAsBool(parameters, "VERBOSE", context)

        if source_in is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        # Transform the in source into a vector layer
        vector_layer_in = source_in.materialize(QgsFeatureRequest(), feedback)

        # Normalize and extract QGS input features
        qgs_features_in, geom_type = Simplify.normalize_in_vector_layer(vector_layer_in, feedback)

        # Validate input geometry type
        if geom_type not in (QgsWkbTypes.LineString, QgsWkbTypes.Polygon):
            raise QgsProcessingException("Can only process: (Multi)LineString or (Multi)Polygon vector layers")

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                               vector_layer_in.fields(),
                                               geom_type,
                                               vector_layer_in.sourceCrs())

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))

        # Set progress bar to 1%
        feedback.setProgress(1)

        # Call ReduceBend algorithm
        rb_return = Simplify.reduce(qgs_features_in, tolerance, validate_structure, feedback)

        for qgs_feature_out in rb_return.qgs_features_out:
            sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

        # Push some output statistics
        feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
        feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
        feedback.pushInfo("Number of iteration needed: {0}".format(rb_return.nbr_pass))
        feedback.pushInfo("Number of vertice deleted: {0}".format(rb_return.nbr_vertice_deleted))
        if validate_structure:
            if rb_return.is_structure_valid:
                status = "Valid"
            else:
                status = "Invalid"
            feedback.pushInfo("Debug - State of the internal data structure: {0}".format(status))
        if verbose:
            for i in range(rb_return.nbr_pass):
                str_value = "Iteration: {}; Bends detected: {}; Bend reduced: {}" \
                            .format(i, rb_return.nbr_bend_detected[i], rb_return.nbr_bend_reduced[i])
                feedback.pushInfo("Verbose - {0}".format(str_value))

        return {"OUTPUT": dest_id}


# --------------------------------------------------------
# Start of the algorithm
# --------------------------------------------------------

# Define global constant


class RbResults:
    """Class defining the stats and result"""

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'nbr_vertice_deleted',  'qgs_features_out', 'nbr_pass',
                 'is_structure_valid')

    def __init__(self):
        """Constructor that initialize a RbResult object.

        :param: None
        :return: None
        :rtype: None
        """

        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_vertice_deleted = 0
        self.qgs_features_out = None
        self.nbr_pass = 0
        self.is_structure_valid = None


class Simplify:
    """Main class for bend reduction"""

    @staticmethod
    def normalize_in_vector_layer(in_vector_layer, feedback):
        """Method used to normalize the input vector layer

        Two processing are used to normalized the input vector layer
         - execute "Multi to single part" processing in order to accept even multipolygon
         - execute "Drop  Z and M values" processing as they are not useful
         - Validate if the resulting layer is Point LineString or Polygon

        :param in_vector_layer:  Input vector layer to normalize
        :param feedback: QgsFeedback handle used to communicate with QGIS
        :return Output vector layer and Output geometry type
        :rtype Tuple of 2 values
        """

        # Execute MultiToSinglePart processing
        feedback.pushInfo("Start normalizing input layer")
        params = {'INPUT': in_vector_layer,
                  'OUTPUT': 'memory:'}
        result_ms = processing.run("native:multiparttosingleparts", params, feedback=feedback)
        ms_part_layer = result_ms['OUTPUT']

        # Execute Drop Z M processing
        params = {'INPUT': ms_part_layer,
                  'DROP_M_VALUES': True,
                  'DROP_Z_VALUES': True,
                  'OUTPUT': 'memory:'}
        result_drop_zm = processing.run("native:dropmzvalues", params, feedback=feedback)
        drop_zm_layer = result_drop_zm['OUTPUT']

        # Extract the QgsFeature from the vector layer
        qgs_in_features = []
        qgs_features = drop_zm_layer.getFeatures()
        for qgs_feature in qgs_features:
            qgs_in_features.append(qgs_feature)
        if len(qgs_in_features) > 1:
            geom_type = qgs_in_features[0].geometry().wkbType()
        else:
            geom_type = drop_zm_layer.wkbType()  # In case of empty layer
        feedback.pushInfo("End normalizing input layer")

        return qgs_in_features, geom_type

    @staticmethod
    def douglas_peucker(qgs_in_features, tolerance, validate_structure=False, feedback=None):
        """Main static method used to launch the bend reduction.

        :param: [QgsFeatures] qgs_features: List of features to process.
        :param: tolerance: Simplification tolerance.
        :param: Bool validate_structure: Validate internal data structure after processing (for debugging)
        :param: QgsFeedback feedback: Handle for interaction with QGIS.
        :return: Statistics and result object.
        :rtype: RbResult
        """

        dp = Simplify(qgs_in_features, tolerance, validate_structure, feedback)
        results = dp.reduce()

        return results

    @staticmethod
    def create_polygon(i, j, qgs_points):
        """This method create a polygon from a subset of point (vertice)

        Note: The first vertice is also the last vertice

        :param: i: Start of point in the list
        :param: j: End of point in the list
        :param: qgs_points: List of QgsPoint
        :return: A polygon formed by a subset of the list of QgsPoint
        :rtype: QgsPolygon
        """

        # Create the list of point to create
        if i < j:
            index = list(range(i, j+1)) + [i]
        else:
            index = list(range(i, len(qgs_points))) + list(range(0, j+1)) + [i]  # Manage circular array

        qgs_sub_points = [qgs_points[k] for k in index]
        qgs_polygon = QgsPolygon(QgsLineString(qgs_sub_points))

        return qgs_polygon

    @staticmethod
    def validate_simplicity(qgs_geoms_with_itself, qgs_geom_new_subline):
        """Validate the simplicitity constraint

        This constraint assure that the new sub line is not intersecting with any other segment of the same line

        :param: qgs_geoms_with_itself: List of QgsLineString segment to verify for self intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        geom_engine_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet().clone())
        for qgs_geom_potential in qgs_geoms_with_itself:
            de_9im_pattern = geom_engine_subline.relate(qgs_geom_potential.constGet().clone())
            # de_9im_pattern[0] == '0' means that their interiors intersect (crosses)
            # de_9im_pattern[1] == '0' means that one extremity is touching the interior of the other (touches)
            if de_9im_pattern[0] == '0' or de_9im_pattern[1] == '0':
                # The new sub line intersect or touch with itself. The result would create a non OGC simple line
                constraints_valid = False
                break

        return constraints_valid

    @staticmethod
    def validate_intersection(qgs_geom_with_others, qgs_geom_new_subline):
        """Validate the intersection constraint

        This constraint assure that the new sub line is not intersecting with any other lines (not itself)

        :param: qgs_geoms_with_others: List of QgsLineString segment to verify for intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        for qgs_geom_potential in qgs_geom_with_others:
            if not qgs_geom_potential.disjoint(qgs_geom_new_subline):
                # The bend area intersects with a point
                constraints_valid = False
                break

        return constraints_valid

    @staticmethod
    def validate_sidedness(qgs_geom_with_others, qgs_geom_bend):
        """Validate the sidedness constraint

        This constraint assure that the new sub line will not change the relative position of an object compared to
        the polygon formed by the bend to reduce. ex.: an interior ring of a polygon going outside of the exterior ring.

        :param: qgs_geoms_with_others: List of QgsLineString segment to verify for intersection
        :param: qgs_geom_bend: QgsPolygon formed by the bend to reduce
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        for qgs_geom_potential in qgs_geom_with_others:
            if qgs_geom_bend.contains(qgs_geom_potential):
                # A feature is totally located inside
                constraints_valid = False
                break

        return constraints_valid

    __slots__ = ('qgs_in_features', 'tolerance', 'validate_structure', 'feedback', 'rb_collection', 'eps',
                 'rb_results', 'rb_features', 'rb_geoms',
                 'bends_reduced')

    def __init__(self, qgs_in_features, tolerance, validate_structure, feedback):
        """Constructor for the bend reduction.

       :param: qgs_in_features: List of features to process.
       :param: tolerance: Float tolerance of the diameter of the bend to reduce.
       :param: validate_structure: flag to validate internal data structure after processing (for debugging)
       :param: feedback: QgsFeedback handle for interaction with QGIS.
       """

        self.qgs_in_features = qgs_in_features
        self.tolerance = tolerance
        self.validate_structure = validate_structure
        self.feedback = feedback
        self.bends_reduced = []  # List containing the reduced bend
        self.eps = None
        self.rb_results = None
        self.rb_features = None
        self.rb_geoms = None
        self.rb_collection = None

    def reduce(self):
        """Main method to reduce line string.

        :return: Statistics and result object.
        :rtype: RbResult
        """

        """
        #  Code used for the profiler (uncomment if needed)
        import cProfile, pstats, io
        from pstats import SortKey
        pr = cProfile.Profile()
        pr.enable()
        """

        # Calculates the epsilon and initialize some stats and results value
        self.eps = Epsilon(self.qgs_in_features)
        self.eps.set_class_variables()
        self.rb_results = RbResults()

        # Create the list of GsPolygon, GsLineString and GsPoint to process
        self.rb_features = self.create_rb_feature()
        self.rb_results.in_nbr_features = len(self.qgs_in_features)

        # Pre process the LineString: remove to close point and co-linear points
        self.rb_geoms = self.pre_simplification_process()

        # Create the GsCollection a spatial index to accelerate search
        self.rb_collection = GsCollection()
        self.rb_collection.add_features(self.rb_geoms)

        # Execute the line simplification for each LineString
        self._simplify_lines()

        # Recreate the QgsFeature
        qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in self.rb_features]

        # Set return values
        self.rb_results.out_nbr_features = len(qgs_features_out)
        self.rb_results.qgs_features_out = qgs_features_out

        # Validate inner spatial structure. For debug purpose only
        if self.rb_results.is_structure_valid:
            self.rb_collection.validate_integrity(self.rb_geoms)

        #  Code used for the profiler (uncomment if needed)
        """
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        """

        return self.rb_results

    def create_rb_feature(self):
        """Create the different GsFeatures from the QgsFeatures.

        :return: List of rb_features
        :rtype: [GsFeature]
        """

        rb_features = []

        for qgs_feature in self.qgs_in_features:
            qgs_geom = qgs_feature.geometry()  # extract the Geometry

            if GsFeature.is_polygon(qgs_geom.wkbType()):
                rb_features.append(GsPolygon(qgs_feature))
            elif GsFeature.is_line_string(qgs_geom.wkbType()):
                rb_features.append(GsLineString(qgs_feature))
            elif GsFeature.is_point(qgs_geom.wkbType()):
                rb_features.append(GsPoint(qgs_feature))
            else:
                raise QgsProcessingException("Internal geometry error")

        return rb_features

    def pre_simplification_process(self):
        """This method execute the pre simplification process

        Pre simplification process applies only to closed line string and is used to find the 2 points that are
        the distant from each other using the oriented bounding box

        :return: List of rb_geom
        :rtype: [RbGeom]
        """

        # Create the list of RbGeom ==> List of geometry to simplify
        sim_geoms = []
        for rb_feature in self.rb_features:
            sim_geoms += rb_feature.get_rb_geom()

        # Find the 2 most distant vertice in a closed geometry
        for sim_geom in sim_geoms:
            qgs_line_string = sim_geom.qgs_geom.constGet()
            if qgs_line_string.isClosed():
                qgs_oriented_bbox = qgs_line_string.orientedMinimumBoundingBox()
                qgs_geom_bbox = qgs_oriented_bbox[0]
                qgs_points = qgs_geom_bbox.constGet().exteriorRing().points()  # Extract vertice of the bounding box
                length_axis_0 = qgs_points[0].distance(qgs_points[1])
                length_axis_1 = qgs_points[1].distance(qgs_points[2])
                if length_axis_0 < length_axis_1:
                    # Find the two shortest sides of the bounding box
                    line0 = QgsLineString(qgs_points[0], qgs_points[3])
                    line1 = QgsLineString(qgs_points[0], qgs_points[1])
                    qgs_points = qgs_line_string.points()
                    distances_line0 = [line0.distance(qgs_point) for qgs_point in qgs_points]
                    distances_line1 = [line1.distance(qgs_point) for qgs_point in qgs_points]
                    new_start_end = distances_line0.index(max(distances_line0))
                    sim_geom.furthest_index = distances_line1.index(max(distances_line1))
                    if new_start_end == sim_geom.furthest_index:
                        # Special case should not happen but just in case
                        new_start_end = 0
                        sim_geom.furthest_index = len(qgs_points)//2
                if new_start_end != 0:
                    # Move the first/last vertice to the new location
                    new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end + 1]
                    sim_geom.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))

        return sim_geoms

    def _simplify_lines(self):
        """Loop over the geometry until there is no more bend to reduce

        An iterative process for bend reduction is needed in order to maximise the bend reduction.  The process
        will always stabilize and exit when there are no more bends to reduce.

        """

        min_nbr_pass = 2
        nbr_geoms = 100.0 / len(self.rb_geoms) if len(self.rb_geoms) >= 1 else 0
        while True:
            self.feedback.setProgress(max(1, int(self.count_rb_geoms_done() * nbr_geoms)))
            nbr_vertice_deleted = 0
#            for rb_geom in self.rb_geoms:
#                if self.feedback.isCanceled():
#                    break
#                if not rb_geom.is_simplest:  # Only process geometry that are not at simplest form
#                    self.delete_co_linear(rb_geom)
#                    nbr_bend_detected = ReduceBend.detect_bends(rb_geom)
#                    if rb_geom.need_pivot:
#                        #  Pivoting a closed line moves the first/last vertice on a bend that do not need simplification
#                        ReduceBend.pivot_closed_line(rb_geom, self.diameter_tol)
#                        nbr_bend_detected = ReduceBend.detect_bends(rb_geom)  # Bend detection needed after pivot
#                    ReduceBend.flag_bend_to_reduce(rb_geom, self.diameter_tol)
#                    nbr_bend_reduced += self.process_bends(rb_geom)

            self.rb_results.nbr_bend_reduced.append(nbr_vertice_deleted)

            # While loop breaking condition
            if self.rb_results.nbr_pass > min_nbr_pass and nbr_vertice_deleted == 0:
                break
            self.rb_results.nbr_pass += 1

        return

    def count_rb_geoms_done(self):
        """Count the number of geometry  that are at there simplest form

        """

        nbr_done = 0
        for rb_geom in self.rb_geoms:
            if rb_geom.is_simplest:
                nbr_done += 1

        return nbr_done

#    def delete_co_linear(self, rb_geom):
#        """Delete co-linear vertice on a LineString
#
#        This method delete co-linear and near co-linear vertice because they are unnecessary but moreover in certain
#        condition when they are forming near 0 (empty) area these vertices are creating spatial calculus errors
#
#        :param: RbGeom rb_geom: Geometry to delete co-linear vertices
#        """
#
#        # Build the list of angles for each vertice
#        vertex_ids_to_del = []
#        angles = ReduceBend.get_angles(rb_geom.qgs_geom.constGet())
#        if rb_geom.qgs_geom.constGet().isClosed() and len(angles) >= 1:
#            del angles[0]  # Do not process the start/end vertice (even if co-linear)
#        for i, angle in enumerate(angles):
#            if abs(angle - math.pi) <= Epsilon.ZERO_ANGLE or abs(angle) <= Epsilon.ZERO_ANGLE:
#                # Co-linear point or flat angle delete the current point
#                vertex_ids_to_del.append(i+1)
#
#        # Delete co-linear vertex
#        for vertex_id_to_del in reversed(vertex_ids_to_del):
#            self.rb_collection.delete_vertex(rb_geom, vertex_id_to_del, vertex_id_to_del)
#
#        # Special case to process closed line string to find ans delete co-linear points at the first/last vertice
#        if rb_geom.qgs_geom.constGet().isClosed():
#            num_points = rb_geom.qgs_geom.constGet().numPoints()
#            if num_points >= 5:  # Minimum of 5 vertices are needed to have co-linear vertices in closed line
#                qgs_ls = QgsLineString([rb_geom.qgs_geom.vertexAt(num_points-2),
#                                        rb_geom.qgs_geom.vertexAt(0),
#                                        rb_geom.qgs_geom.vertexAt(1)])
#                angles = ReduceBend.get_angles(qgs_ls)
#                angle = angles[0]
#                if abs(angle - math.pi) <= Epsilon.ZERO_ANGLE or abs(angle) <= Epsilon.ZERO_ANGLE:
#                    self.rb_collection.delete_vertex(rb_geom, 0, 0)
#
#        if rb_geom.qgs_geom.length() <= Epsilon.ZERO_RELATIVE:
#            # Something wrong.  do not try to simplify the LineString
#            rb_geom.is_simplest = True
#
#        return

    def validate_constraints(self, ind, rb_geom):
        """Validate the spatial relationship in order maintain topological structure

        Three distinct spatial relation are tested in order to assure that each bend reduce will continue to maintain
        the topological structure in a feature between the features:
         - Simplicity: Adequate validation is done to make sure that the bend reduction will not cause the feature
                       to cross  itself.
         - Intersection : Adequate validation is done to make sure that a line from other features will not intersect
                          the bend being reduced
         - Sidedness: Adequate validation is done to make sure that a line is not completely contained in the bend.
                      This situation can happen when a ring in a polygon complete;y lie in a bend ans after bend
                      reduction, the the ring falls outside the polygon which make it invalid.

        Note if the topological structure is wrong before the bend correction no correction will be done on these
        errors.

        :param: ind: Index number of the bend to process
        :param: rb_geom: Geometry used to validate constraints
        :param: detect_alternate_bend: Indicates if alternate bend can be find when self intersection is detected
        :return: Flag indicating if the spatial constraints are valid for this bend reduction
        :rtype: Bool
        """

        constraints_valid = True
        bend = rb_geom.bends[ind]
        b_box = bend.qgs_geom_bend.boundingBox()
        qgs_geoms_with_itself, qgs_geoms_with_others = \
            self.rb_collection.get_segment_intersect(rb_geom.id, b_box, bend.qgs_geom_old_subline)

        # First: check if the bend reduce line string is an OGC simple line
        # We test with a tiny smaller line to ease the testing and false positive error
        if bend.qgs_geom_new_subline.length() >= Epsilon.ZERO_RELATIVE:
            constraints_valid = Simplify.validate_simplicity(qgs_geoms_with_itself, bend.qgs_geom_new_subline)
        else:
            # Error in the input file
            qgs_line_string = bend.qgs_geom_new_subline.constGet()
            x = qgs_line_string.startPoint().x()
            y = qgs_line_string.startPoint().y()
            text = "Possibly non OGC simple feature at {},{} use Fix Geometries".format(x, y)
            self.feedback.pushInfo(text)

        # Second: check that the new line does not intersect any other line or points
        if constraints_valid:
            constraints_valid = Simplify.validate_intersection(qgs_geoms_with_others, bend.qgs_geom_new_subline)

        # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
        # sidedness or relative position error
        if constraints_valid:
            constraints_valid = Simplify.validate_sidedness(qgs_geoms_with_others, bend.qgs_geom_bend)

        return constraints_valid

    def process_line(self, line, pass_nbr):
        """
        This method is simplifying a line with the Douglas Peucker algorithm plus contraints checking if they are enabled.

        This method is checking the line differently the first time from the remaining time.  The idea behind it is only to
        have a faster process. We assume that most of the lines will not have a problem so we check for problems (SIMPLE_LINE,
        CROSSING_LINE and SIDEDNESS) against the whole line for the first time. If there are some problems tthe next time we will
        check each sub portion of the line. This strategy is making a huge difference in time.

        Parameters:
            line: The line to process
            pass_nbr: The number of the pass. At their first pass we process the line as a whole. We do not
            check each sub portion of line simplification. For the other passes we check each sub portion of the line
            for constraints

        Return value:
            True: The line is simplified
            False: The line is not simplified
        """

#        index = set()  # Set of the points to keep
#        stack = []  # Stack to simulate the recursion
#        line.simpliest = True
#        first = 0
#        last = len(line.coords_dual) - 1
#        stack.append((first, last))

#        while stack:
#            (first, last) = stack.pop()
#            if first + 1 < last:  # The segment to check has only 2 points
#                add_point = True
#                (farthest_index, farthest_dist) = self._find_farthest_point(line, first, last)
#                if farthest_index != first:
#                    if farthest_dist <= line.ma_properties[_TOLERANCE]:

#                        if (pass_nbr != 0):
#                            # For all the pass except the first one we check each sub line simplification individually
#                            line_simple_line = LineString(line._coords_dual[:first + 1] + line.coords_dual[last:])
#                            new_segment_coords = [line.coords_dual[first], line.coords_dual[last]]
#                            old_segment_coords = line.coords_dual[first:last + 1]
#                            line_crossing_line = LineString(new_segment_coords)
#                            sidedness_polygon = GenUtil.calculate_sidedness_polygon(LineString(old_segment_coords),
#                                                                                    LineString(new_segment_coords))
#
#                            conflict_type = GenUtil.test_constraints(self, None, line_simple_line, line_crossing_line,
#                                                                     sidedness_polygon, self.s_container, line._sci_id)
#                        else:
#                            # We check for conflict only at the end of the process so here we assume no conflict
#                            conflict_type = None
#
#                        if (conflict_type is not None):
#                            line.simpliest = False  # This line is not at it's simplest form since
#                        else:  # a constraint is blocking the simplification
#                            index.update([first, last])
#                            add_point = False
#
#                    if add_point:
#                        stack.append((first, farthest_index))
#                        stack.append((farthest_index, last))
#                else:
#                    index.update([first, last])
#
#            else:
#                index.update([first, last])
#
#        replacement_index = list(index)
#
#        if (line.is_closed and (len(replacement_index) <= 3)):
#            #           Closed line must have at least 4 vertices
#            replacement_index = self._process_closed_line(line)
#
#        # Check if the line has been simplified
#        nbr_vertice_simplified = len(line.coords_dual) - len(replacement_index)
#        if nbr_vertice_simplified == 0:
#            simplified = False  # No change done (same quantity of coordinates)
#            line.is_simplest = True  # The line is at its simplest form
#        else:
#            new_coords = [line.coords_dual[i] for i in sorted(replacement_index)]
#            if (pass_nbr != 0):
#                # If we process each sub modifification of the line inividually
#                simplified = True  # The line has been simplified
#            else:
#                # For the first iteration we process the line as a whole
#                # Check for conglict
#                line_simple_line = LineString(new_coords)
#                new_segment_coords = new_coords
#                old_segment_coords = line.coords_dual
#                line_crossing_line = LineString(new_segment_coords)
#                sidedness_polygon = GenUtil.calculate_sidedness_polygon(LineString(old_segment_coords),
#                                                                        LineString(new_segment_coords))
#
#                conflict_type = GenUtil.test_constraints(self, None, line_simple_line, line_crossing_line,
#                                                         sidedness_polygon, self.s_container, line._sci_id)
#                if (conflict_type is None):
#                    simplified = True  # The line was  simplified
#                    line.is_simplest = True  # If at the first pass the whole line as no conflict it is at its simplest form
#                else:
#                    simplified = False  # The line was not simplified
#
#            if (simplified):
#                for i in xrange(nbr_vertice_simplified): self.stats.add_stats(_ALGO)
#                line.update_coords(new_coords, self.s_container)
#
#        return simplified