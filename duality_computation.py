import numpy as np
from matplotlib import pyplot as plt
import random
import math

########################################
#### Methods for duality transforms ####
########################################

'''
Inversion duality transform
Input: A circle through the origin defined by a center point (x, y)
Returns: The dual line defined by a tuple (slope, intercept, orientation)
'''
def duality1_circleToLine(circle_center):
    x, y = circle_center
    slope = -x/y
    intercept = 1/(2*y)
    orientation = 1 if y > 0 else -1
    return (slope, intercept, orientation)

'''
Inversion duality transform: translation back to primal
Input: Two points in the form (x, y) that define a line segment
Returns: The arc (x, y, r, theta0, theta1) in the primal
         theta0, theta1 are measured in degrees
'''
def inverseDuality1_segmentToArc(point1, point2):
    (x0, y0), (x1, y1) = point1, point2

    # Find primal circle by first converting to point-slope form
    slope = (y1 - y0)/(x1 - x0)
    intercept = y0 - slope * x0
    circle_y = 1/(2*intercept)
    circle_x = -slope * circle_y
    circle_radius = math.sqrt(circle_x ** 2 + circle_y ** 2)

    # Find bounding angles
    primal_x0, primal_y0 = x0/(x0 ** 2 + y0 ** 2), y0/(x0 ** 2 + y0 ** 2)
    primal_x1, primal_y1 = x1/(x1 ** 2 + y1 ** 2), y1/(x1 ** 2 + y1 ** 2)

    theta0 = math.atan((primal_y0 - circle_y)/(primal_x0 - circle_x)) * 180/math.pi
    theta1 = math.atan((primal_y1 - circle_y)/(primal_x1 - circle_x)) * 180/math.pi

    # Deal with the fact that atan returns a number in [-pi/2, pi/2], which cannot differentiate between quadrants
    if   primal_x0 - circle_x < 0: theta0 = 180 + theta0
    if   primal_x1 - circle_x < 0: theta1 = 180 + theta1

    return (circle_x, circle_y, circle_radius, theta0, theta1)

'''
Point-line dual transformation
Input: tuple (m, b) representing the line y = mx + b
Output: tuple (x, y) representing a point
'''
def duality2_lineToPoint(line):
    m, b, orientation = line
    return (-m, b)

'''
Point-line dual transformation: translation back to primal
Input: tuple (x, y) representing a point
Output: tuple (m, b) representing the line y = mx + b
'''
def inverseDuality2_pointToLine(point):
    x, y = point
    return (-x, y)

#################################
### Solve convex hull problem ###
#################################

'''
convex_hull_step
    Processes one more point with the sweepline convex hull algorithm
Input:
- points: list of points (x, y) sorted by x-coordinate
- hull_type: 'upper' or 'lower'
- state: tuple (hull, step_num, done) where
         hull: list containing points (x, y) in the upper/lower hull (depending on 'hull_type')
         step_num: int, the index of the next point to process
         done: boolean that is True if and only if all points have been processed
Output:
- state: tuple like the input parameter, representing the state after the method terminates
'''
def convex_hull_step(points, hull_type, state = None):
    if state is None:
        if len(points) == 0: return [], [], None, True
        hull = [points[0]]
        next_step = 1
        done = (next_step >= len(points))
        state = (hull, next_step, done)
        return state

    hull, step_num, done = state
    curr_point = points[step_num]

    if hull_type == 'upper':
        while len(hull) >= 2 and left_turn(hull[-2], hull[-1], curr_point):
            hull = hull[:-1]    
        hull.append(curr_point)
    elif hull_type == 'lower':
        while len(hull) >= 2 and not left_turn(hull[-2], hull[-1], curr_point):
            hull = hull[:-1]
        hull.append(curr_point)

    next_step = step_num + 1
    done = (next_step >= len(points))
    new_state = (hull, next_step, done)
    return new_state


def left_turn(p1, p2, p3):
    determinant = np.linalg.det([[p1[0], p1[1], 1], [p2[0], p2[1], 1], [p3[0], p3[1], 1]])
    return determinant > 0

#########################################################
### Translate into language of halfplane intersection ###
#########################################################

'''
Find the interpretation of the current convex hull in space of halfplanes
Returns a list of points (x, y) that define the intersection

Note:
  If hull_type is 'lower', the input 'convex_hull_points' must be sorted by increasing x-coordinate
                           and the method returns points in order of increasing x-coordinate
                           (ie. clockwise order around the region)
  If hull_type is 'upper', points are sorted by decreasing x-coordinate
                           (ie. clockwise order around the region)
'''
def halfplane_envelope(convex_hull_points, type):
    lines = [inverseDuality2_pointToLine(point) for point in convex_hull_points]
    intersections = neighboring_intersections(lines)
    
     # Set min_x, max_x as proxies for +/- infinity
     # Differentiate upper/lower hull so that no points have the same x-coordinate
    if type == 'upper':
        min_x, max_x = -1000, 1000
    else:
        min_x, max_x = -1010, 1010
    
    if type == 'upper':
        m_right, b_right = lines[0]
        m_left, b_left = lines[-1]
        right_y = m_right * max_x + b_right
        left_y = m_left * min_x + b_left
        intersections.insert(0, (max_x, right_y))
        intersections.append((min_x, left_y))
    elif type == 'lower':
        m_right, b_right = lines[-1]
        m_left, b_left = lines[0]
        right_y = m_right * max_x + b_right
        left_y = m_left * min_x + b_left
        intersections.insert(0, (min_x, left_y))
        intersections.append((max_x, right_y))

    return intersections

'''
Merge two convex chains that are guaranteed to intersect at at most 1 point
Returns a list of intersection points (x, y) that traverse the intersection region
    in order of descending y-coordinate
'''
def merge_halfplanes(upper_envelope, lower_envelope):
    if upper_envelope is None or len(upper_envelope) == 0:
        edges = lower_envelope
    elif lower_envelope is None or len(lower_envelope) == 0:
        edges = upper_envelope
    else:
        upper_envelope.reverse()
        intersection, (lower_index, upper_index) = intersection_of_envelopes(upper_envelope, lower_envelope)
        
        if intersection is None: return []
        intersection_left = (upper_envelope[0][1] < lower_envelope[0][1])

        if intersection_left:
            lower_edges = lower_envelope[:lower_index + 1]
            upper_edges = upper_envelope[:upper_index + 1]
            upper_edges.reverse()

            edges = lower_edges
            edges.append(intersection)
            edges.extend(upper_edges)
        else:
            lower_edges = lower_envelope[lower_index + 1:]
            upper_edges = upper_envelope[upper_index + 1:]
            upper_edges.reverse()

            edges = upper_edges
            edges.append(intersection)
            edges.extend(lower_edges)
        
    return edges

## Helper methods for halfplane methods ##

'''
Finds (at most one) intersection between two convex chains
Input: two chains of points (x, y) in order of increasing x-coordinate
Returns:
  If an intersection is found, return (x, y), (upper_index, lower_indes)
    where (x, y) is the location of the intersection
    and upper_index, lower_index are the highest indices of points *before* the intersection
  If no intersection is found, return None, (None, None)
'''
def intersection_of_envelopes(upper_envelope, lower_envelope):
    stopping_points = [(x, y, 'upper') for x, y in upper_envelope[1:]]
    stopping_points.extend([(x, y, 'lower') for x, y in lower_envelope[1:]])
    stopping_points.sort(key = lambda point: point[0])

    envelopes = [lower_envelope, upper_envelope]
    curr_index = [0, 0]

    # Initial check for intersection
    lower_p0, lower_p1 = lower_envelope[0], lower_envelope[1]
    upper_p0, upper_p1 = upper_envelope[0], upper_envelope[1]
    intersection = intersection_of_segments(lower_p0, lower_p1, upper_p0, upper_p1)
    if intersection is not None:
        return intersection, curr_index

    for i, point in enumerate(stopping_points):
        x, y, point_type = point
        if point_type == 'intersection':
            return (x, y), curr_index

        primary_index = 1 if point_type == 'upper' else 0
        secondary_index = 0 if point_type == 'upper' else 1
        primary_label = 'upper' if primary_index == 1 else 'lower'
        secondary_label = 'lower' if primary_index == 1 else 'upper'
        primary_envelope = envelopes[primary_index]
        secondary_envelope = envelopes[secondary_index]
        currently_testing = curr_index[primary_index] + 1

        # If we need to process the last points (ie at infinity),
        # we know there can be no intersection
        if (currently_testing + 1 >= len(primary_envelope)
            or curr_index[secondary_index] + 1 >= len(secondary_envelope)):
            return None, [None, None]

        assert((x,y) == primary_envelope[currently_testing])
        primary_point1, primary_point2 = ((x, y), primary_envelope[currently_testing + 1])
        secondary_point1, secondary_point2 = (secondary_envelope[curr_index[secondary_index]], secondary_envelope[curr_index[secondary_index] + 1])
        intersection = intersection_of_segments(primary_point1, primary_point2, secondary_point1, secondary_point2)

        curr_index[primary_index] += 1

        if intersection is not None:
            stopping_points.append((intersection[0], intersection[1], 'intersection'))
            stopping_points.sort()

'''
Returns a list (x, y) of intersection points between a lines next to each other
in a list of (slope, intercept)
'''
def neighboring_intersections(lines):
    x_list = []
    y_list = []
    for i in range(len(lines) - 1):
        m_1, b_1 = lines[i]
        m_2, b_2 = lines[i + 1]
        intersection_x = (b_2 - b_1)/(m_1 - m_2)
        intersection_y = m_1 * intersection_x + b_1
        x_list.append(intersection_x)
        y_list.append(intersection_y)
    intersections = [(x, y) for x, y in zip(x_list, y_list)]
    return intersections

'''
Finds the slope and intercept of a line through two points (x, y)
'''
def convert_to_slope_intercept(point1, point2):
    x0, y0 = point1
    x1, y1 = point2
    m = (y1 - y0)/(x1 - x0)
    b = y0 - m * x0
    return (m,b)

'''
Finds the intersection point between to line segments, returning None if none exist
'''
def intersection_of_segments(line1_point1, line1_point2, line2_point1, line2_point2):
    line1 = convert_to_slope_intercept(line1_point1, line1_point2)
    line2 = convert_to_slope_intercept(line2_point1, line2_point2)
    possible_intersection = intersection_points([line1, line2])
    if len(possible_intersection) == 0: return None
    
    x, y = possible_intersection[0]
    if (line1_point1[0] < x and x < line1_point2[0] and
        line2_point1[0] < x and x < line2_point2[0]):
        return (x,y)
    else:
        return None

'''
Returns a list (x, y) of intersection points between a list of lines (slope, intercept)
'''
def intersection_points(lines):
    x_list = []
    y_list = []
    for i in range(len(lines)):
        for j in range(len(lines)):
            if i == j: continue
            m_1, b_1 = lines[i]
            m_2, b_2 = lines[j]
            intersection_x = (b_2 - b_1)/(m_1 - m_2)
            intersection_y = m_1 * intersection_x + b_1
            x_list.append(intersection_x)
            y_list.append(intersection_y)
    intersections = [(x, y) for x, y in zip(x_list, y_list)]
    return intersections

#####################################################
### Translate into language of  disk intersection ###
#####################################################

'''
Input: a list of (x,y) points defining a halfplane intersection
Returns: a list of (x, y, r, theta0, theta1) arcs
'''
def disk_intersection(halfplane_intersection):
    x_points, y_points = zip(*halfplane_intersection)
    segments = [((x_points[i], y_points[i]), (x_points[i+1], y_points[i+1])) for i in range(len(x_points) - 1)]
    arcs = [inverseDuality1_segmentToArc(point1, point2) for point1, point2 in segments]
    return arcs