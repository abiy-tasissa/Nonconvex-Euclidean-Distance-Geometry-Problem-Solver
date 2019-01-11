# Nonconvex-Euclidean-Distance-Geometry-Problem-Solver
A fast nonconvex algorithm for the Euclidean Distance Geometry Problem. Formally, consider a set of n points where partial inter-point distance information is provided. The goal of Euclidean Distance Geometry problem is to find the coordinate of the points given this partial information. The partial information can be exact or noisy and the solver handles these different cases. These set of codes were written as part of a research project on the Euclidean Distance Geometry Problem. The codes are written by Professor [Rongjie Lai](http://homepages.rpi.edu/~lair/) and Abiy Tasissa. 

## MATLAB files description
`demo_globalrecon.m`: This is the main file which loads different kinds of data file from the data folder and reads the (x,y,z) coordinates of the data. Using these coordinates, the full distance matrix for the data is constructed. Then, some entries of the distance matrix are selected uniformly at random. Once these information and some parameters for the main algorithm are set, the script calls the main algorithms.  

`alternating_completion.m`: This is the main algorithm for the Euclidean Distance Geometry Problem for the case of exact partial information. 

`alternating_completion_noisy.m`: This is the main algorithm for the Euclidean Distance Geometry Problem for the case of noisy partial information. The noise is assumed to be additive Gaussian.

`BBGradient.m`: This implements the BB gradient method coupled with nonmonotone line search. 

`ReadOFF.m, ViewMesh.m, ViewPC.m`: These scripts are used to read different types of the data file and generate points, mesh. They are also used in viewing the outputs of the main algorithms.

## List of data files tested
* `1k.off`: Sphere with 1002 points
* `cow.off`: Cow with 2601 points
* `ptswiss.mat`: Swiss roll data with 2048 points
* `UScities.mat`: Data of 2920 US cities

The first and second data sets are taken from [here](http://visionair.ge.imati.cnr.it/ontologies/shapes/search.jsp).
The third data was obtained by simply plotting the parametric equations of a Swiss roll. The last data uses Latitude
and Longitude of US cities, in different zip codes, to generate the point coordinates. 

## Instructions

The starting script is `demo_globalrecon.m`. Choose a data, sampling rate, set algorithm parameters and call either the exact or noisy solver. 
## References

Abiy Tasissa and Rongjie Lai, "Exact Reconstruction of Euclidean Distance Geometry Problem Using Low-rank Matrix Completion," in IEEE Transactions on Information Theory, 2018. doi: 10.1109/TIT.2018.2881749
[URL](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8537996&isnumber=4667673)


## Feedback

Email your feedback to <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

