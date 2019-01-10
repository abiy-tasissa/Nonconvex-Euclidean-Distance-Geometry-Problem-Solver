# Nonconvex-Euclidean-Distance-Geometry-Problem-Solver
A fast nonconvex algorithm for the Euclidean Distance Geometry Problem. Formally, consider a set of n points where partial inter-point distance information is provided. The goal of Euclidean Distance Geometry problem (EDG in short) is to find the coordinate of the points given this partial information. The partial information can be exact or noisy and the solver handles these different cases. These set of codes were written as part of a research project on the Euclidean Distance Geometry Problem. The codes are written by Professor [Rongjie Lai](http://homepages.rpi.edu/~lair/) and Abiy Tasissa. We appreciate Dr. Jia Li for discussions on algorithms and the code at the early stage of the project.

## MATLAB files description
`demo_globalrecon.m`: This is the main file which loads different kinds of data file (some of them from the data folder) and reads the coordinates. Using these coordinates, the full distance matrix for the data is constructed. Then, some entries of the distance matrix are selected uniformly at random. Once these information and some parameters for the main algorithm are set, the script calls the main algorithms.  

`alternating_completion.m`: This is the main algorithm for the Euclidean Distance Geometry Problem for the case of exact partial information. 

`alternating_completion_noisy.m`: This is the main algorithm for the Euclidean Distance Geometry Problem for the case of noisy partial information. The noise is assumed to be additive and Gaussian.

`BBGradient.m`: This implements the BB gradient method coupled with nonmonotone line search.

`ReadOFF.m, READObjShape.m, ViewMesh.m, ViewPC.m`: These scripts are used to read different types of the data file and generate points, mesh. They are also used in viewing the reconstructed 3D geometry from the main algorithms.

## List of data files tested
* `1k.off`: sphere with 1002 points
* `cow_1.obj`: cow with 2601 points
* `kitten_1.obj`: kitten with 2884 points
* `horse.obj`: horse with 19851 points
* `ptswiss.mat`: Swiss roll data with Euclidean distances. It has 2048 points
* `GeodesicDistswiss.mat`: Swiss roll data with geodesic distance. It has 2048 points.
* `UScities.mat`: Data of 2920 US cities

## Instructions

The starting script is `demo_globalrecon.m`. Choose a data, sampling rate, set algorithm parameters and call either the exact or noisy solver. 
## References

Tasissa, Abiy, and Rongjie Lai. "Exact Reconstruction of Euclidean Distance Geometry Problem Using Low-rank Matrix Completion." arXiv preprint arXiv:1804.04310 (2018).

## Feedback

Email your feedback to <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

