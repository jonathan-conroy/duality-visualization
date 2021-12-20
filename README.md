Comp163 Final Project - Intersection of Disks
Jonathan Conroy
Spring 2020

This visualization explores the use of duality transforms to determine the
intersection region of disks that all pass through the origin. Using the
inversion duality transform, the original problem can be transformed into a
problem about the intersection of halfplanes. Using a point-line duality
transform, the halfplane problem can be transformed into a convex hull problem.

This interactive demo shows how finding a convex hull can be interpreted as
finding the intersection of disks.

To use:
    Run 'python visualization.py'
    Click on locations in the left-most plot to add discs
    Press any key to continue (display the dual plots)
    Press any key to continue (find the upper hull)
    Press any key to continue (find the lower hull)
    Press any key to continue (merge envelopes and display final solution)
    Press any key to exit


Files:
    visualization.py:       The driver of the program. Handles visualization.
    duality_computation.py: Provides the methods that apply computational
                            geometry algorithms to solve the problem
    Duality.pptx:           A short set of slides introducing the problem

External Dependencies:
    numpy, matplotlib

References:
Dobkin, D. & Souvaine, Diane. "Computational Geometry -- A User's Guide."
    Advances in Robotic 1: Algorithmic and Geometric Aspects of Robotics,
    1987, pp. 43-93.

Brown, Kevin Q. "Geometric Transforms for Fast Geometric Algorithms."
    Carnegie-Mellon University Department of Computer Science, 1979.
