# Comp163 Final Project - Intersection of Disks
_Jonathan Conroy_  
_Spring 2020_  

This visualization explores the use of duality transforms to determine the
intersection region of disks that all pass through the origin. Using the
inversion duality transform, the original problem can be transformed into a
problem about the intersection of halfplanes [1]. Using a point-line duality
transform, the halfplane problem can be transformed into a convex hull problem [2].

This interactive demo shows how the incremental approach to finding a convex hull can be interpreted as
finding the intersection of disks.

## To use:
- Run 'python visualization.py'  
- Click on locations in the left-most plot to add discs  
- Press any key to continue (display the dual plots)  
- Press any key to continue (find the upper hull)  
- Press any key to continue (find the lower hull)  
- Press any key to continue (merge envelopes and display final solution)  
- Press any key to exit  


## Files:
- _visualization.py_:       The driver of the program. Handles visualization.
- _duality_computation.py_: Provides the methods that apply computational
                                geometry algorithms to solve the problem

## External Dependencies:  
numpy, matplotlib

## References:
[1] Dobkin, D. & Souvaine, Diane. ["Computational Geometry -- A User's Guide."](http://www.cs.tufts.edu/comp/163/notes05/comp_geom__a_users_guide.pdf)
Advances in Robotic 1: Algorithmic and Geometric Aspects of Robotics,
1987, pp. 43-93.

[2] Brown, Kevin Q. ["Geometric Transforms for Fast Geometric Algorithms."](http://reports-archive.adm.cs.cmu.edu/anon/anon/home/ftp/scan/CMU-CS-80-101.pdf)
Carnegie-Mellon University Department of Computer Science, 1979.
  
