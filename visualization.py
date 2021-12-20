import numpy as np
from matplotlib import pyplot as plt
import random
import matplotlib
import math

from duality_computation import *

########################################
######## Set up visualization ##########
########################################

plt.ion()
fig, (circle_ax, halfplane_ax, hull_ax) = plt.subplots(nrows = 1, ncols = 3, figsize = (14,4))
circle_ax.xaxis.set_visible(False)
circle_ax.yaxis.set_visible(False)
halfplane_ax.xaxis.set_visible(False)
halfplane_ax.yaxis.set_visible(False)
hull_ax.xaxis.set_visible(False)
hull_ax.yaxis.set_visible(False)

# Set up circle image
circle_ax.set_xlim(-3, 3)
circle_ax.set_ylim(-3, 3)
circle_ax.scatter([0], [0])
arc_patches = []
arc_plt = circle_ax.add_collection(matplotlib.collections.PatchCollection([]))

# Set up halfplane image
halfplane_ax.set_xlim(-3, 3)
halfplane_ax.set_ylim(-3, 3)
emph_lines_plt, = halfplane_ax.plot([], [], c = 'b', lw = 3, zorder = 3)

# Set up hull image
hull_ax.set_xlim(-3, 3)
hull_ax.set_ylim(-3, 3)

points_plt = hull_ax.scatter([], [])
hull_plt, = hull_ax.plot([], [])
sweepline = hull_ax.axvline(ymax = 0) # Line starts invisible

########################################
#############  Main method  ############
########################################

def main():
    # Draw initial circles
    center_points = get_input_points()

    # Draw initial halfplane lines
    lines = [duality1_circleToLine(circle) for circle in center_points]
    min_x, max_x, min_y, max_y = scale_plot_for_line_intersections(lines)
    x = np.linspace(min_x, max_x, 10)
    for (m, b, orientation) in lines:
        y = m*x + b
        color = 'b' if orientation == 1 else 'r'
        fill_to = max_y if orientation == 1 else min_y

        halfplane_ax.plot(x, y, color)
        halfplane_ax.fill_between(x, y, fill_to, alpha = 0.05, color = color)

    # Calculate upper hull = upper envelope
    upper_points = [duality2_lineToPoint(line) for line in lines if line[-1] == 1]
    min_x, max_x, min_y, max_y = scale_plot_for_points(upper_points)
    draw_points(upper_points)
    upper_envelope = visualize_convex_hull(upper_points, 'upper', min_x, max_x)
    while not plt.waitforbuttonpress(): pass
    draw_lines([])
    draw_points([])

    # Calculate lower hull = lower envelope
    lower_points = [duality2_lineToPoint(line) for line in lines if line[-1] == -1]
    min_x, max_x, min_y, max_y = scale_plot_for_points(lower_points)
    lower_envelope = visualize_convex_hull(lower_points, 'lower', min_x, max_x)
    while not plt.waitforbuttonpress(): pass

    # Merge halfplanes and display final solution
    display_merge_halfplanes(upper_envelope, lower_envelope)
    while not plt.waitforbuttonpress(): pass

####################################
########  Handle input  ############
####################################

input_points = []
def handle_mouse_click(event):
    global input_points
    x, y = event.xdata, event.ydata
    input_points.append((x, y))
    draw_circle(x,y)

def get_input_points():
    global input_points
    cid = fig.canvas.mpl_connect('button_press_event', handle_mouse_click)
    
    # Pause until a key is pressed (ie. the input phase has finished) 
    while not plt.waitforbuttonpress():
        continue
    fig.canvas.mpl_disconnect(cid)
    input_points = list(set(input_points))
    return input_points

###########################################
#####  Main methods for visualization  ####
###########################################

'''
Main method to visualize the sweepline algorithm and the primal interpretation
Input:
- points: list of (x, y)
- hull_typ: 'upper' or 'lower'
- min_x: float, the x-coordinate where the sweepline starts
- max_x: float, the x-coordinate where the sweepline ends
Returns:
- halfplane_intersection: list of (x, y) intersection points defining the envelope
                          in clockwise order -  ie. sorted by increasing x-coordinate
                          for 'lower' and decreasing x-coordinate for 'upper'
'''
def visualize_convex_hull(points, hull_type, min_x, max_x):
    points.sort(key = lambda point: point[0])
    draw_points(points)
    state = None
    halfplane_intersection = None
    display_circle_intersection([])
    emph_lines_plt.set_data([], [])

    sweepline.set_data([min_x, min_x], [0, 1])
    while not plt.waitforbuttonpress(): pass

    num_intervals = 100
    intervals = np.arange(min_x, max_x, (max_x-min_x)/num_intervals)

    # Sweepline animation
    curr_index = 0
    for x in intervals:
        sweepline.set_data([x, x], [0, 1])
        if curr_index < len(points) and x > points[curr_index][0]:
            state = convex_hull_step(points, hull_type, state)
            hull, _, _ = state
            curr_index += 1

            # Update convex hull subplot
            point_sizes = [70 if i < curr_index else 30 for i in range(len(points))]
            points_plt.set_sizes(point_sizes)
            draw_lines(hull)

            # Update halfplane subplot
            halfplane_intersection = display_halfplane_envelope(hull, hull_type)

            # Update circle subplot
            color = 'b' if hull_type == 'upper' else 'r'
            display_circle_intersection(halfplane_intersection, color)
        plt.pause(0.005)
    return halfplane_intersection
    
'''
Displays the interpretation of the convex hull as a halfplane intersection region
Input:
- convex_hull_points: list of (x, y) points
- hull_type: 'upper' or 'lower'
Returns:
- intersection_points: list of (x, y) defining the intersection region in clockwise order
'''
def display_halfplane_envelope(convex_hull_points, hull_type):
    intersection_points = halfplane_envelope(convex_hull_points, hull_type)

    x_points, y_points = zip(*intersection_points)
    color = 'b' if hull_type == 'upper' else 'r'
    emph_lines_plt.set_color(color)
    emph_lines_plt.set_data(x_points, y_points)
    return intersection_points

'''
Displays the interpretation of the intersection of halfplanes as an intersection of disks
Input:
- halfplane_intersection: list of (x, y) defining the intersection region in clockwise order
- color: string representing a matplotlib color
Returns: None
'''
def display_circle_intersection(halfplane_intersection, color = 'b'):
    global arc_patches
    for arc in arc_patches:
        arc.remove()
    arc_patches = []

    if len(halfplane_intersection) == 0: return
    new_arcs = disk_intersection(halfplane_intersection)
    for (x, y, r, theta1, theta2) in new_arcs:
        new_arc = matplotlib.patches.Arc((x, y), 2*r, 2*r, theta1 = theta1, theta2 = theta2, fill = False, ec = color, lw = 3)
        circle_ax.add_patch(new_arc)
        arc_patches.append(new_arc)

'''
Merges halfplanes and displays the final intersection region 
Input:
- upper_envelope: list of (x, y) intersection points defining the envelope in clockwise order
- lower_envelope: list of (x, y) intersection points defining the envelope in clockwise order
'''
def display_merge_halfplanes(upper_envelope, lower_envelope):
    emph_lines_plt.set_data([], [])
    display_circle_intersection([])
    
    edges = merge_halfplanes(upper_envelope, lower_envelope)
    if len(edges) == 0: return
    x, y = zip(*edges)
    halfplane_ax.plot(x, y, c = 'g', lw = 3)

    display_circle_intersection(edges, 'g')



#############################################
#####  Helper methods for visualization  ####
#############################################

'''
Rescale halfplane axis to show all intersections
'''
def scale_plot_for_line_intersections(halfplanes):
    if len(halfplanes) == 1:
        min_x = -1
        max_x = 1
        min_y = -1
        max_y = 1
    else:
        lines = [(m,b) for (m, b, orientation) in halfplanes]
        intersections = intersection_points(lines)
        x_list = [x for (x,y) in intersections]
        y_list = [y for (x,y) in intersections]

        min_x = min(x_list)
        max_x = max(x_list)
        min_y = min(y_list)
        max_y = max(y_list)
        
        x_buffer = 0.2 * abs(max_x - min_x) + 1
        y_buffer = 0.2 * abs(max_y - min_y) + 1
        min_x = min_x - x_buffer
        max_x = max_x + x_buffer
        min_y = min_y - y_buffer
        max_y = max_y + y_buffer

    halfplane_ax.set_xlim(min_x, max_x)
    halfplane_ax.set_ylim(min_y, max_y)
    return min_x, max_x, min_y, max_y

'''
Rescale convex hull axis to show all intersections
'''
def scale_plot_for_points(points):
    if len(points) == 0:
        min_x = -1
        max_x = 1
        min_y = -1
        max_y = 1
    else:
        min_x = min(points, key = lambda point: point[0])[0]
        max_x = max(points, key = lambda point: point[0])[0]
        min_y = min(points, key = lambda point: point[1])[1]
        max_y = max(points, key = lambda point: point[1])[1]

        x_buffer = 0.2 * abs(max_x - min_x) + 1
        y_buffer = 0.2 * abs(max_y - min_y) + 1
        min_x = min_x - x_buffer
        max_x = max_x + x_buffer
        min_y = min_y - y_buffer
        max_y = max_y + y_buffer

    hull_ax.set_xlim(min_x, max_x)
    hull_ax.set_ylim(min_y, max_y)
    return min_x, max_x, min_y, max_y

'''
Display points on convex hull axis
'''
def draw_points(points):
    points_x = [point[0] for point in points]
    points_y = [point[1] for point in points]
    points_plt.set_offsets(np.c_[points_x, points_y])
    points_plt.set_sizes([30 for _ in points])
    plt.draw()
    plt.pause(.1)

'''
Display lines on convex hull axis
The input is a list of points (x, y)
'''
def draw_lines(hull):
    if len(hull) == 0: x, y = [], []
    else: x, y = zip(*hull)
    hull_plt.set_data(x, y)

'''
Display a circle on the convex hull
'''
def draw_circle(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    circle = plt.Circle((x,y), r, fill = False)
    circle_ax.add_artist(circle)


# Call main method
main()
