In this folder, there are some pts files. 
They contain a single branch, from (-0.5, 0, 0) to (0.5, 0, 0).
They are discretized with a number of points equal to the number in the name of the file, +1 (e.g. segment160.pts has 161 points and 160 elements)

BEWARE: in the convergence tests, they are used in smaller 3d domains: (-0.25, 0.25)x(-1,1)^2. Therefore, there is some part of the branch outside the 3d domain (in this way, we have uniformity along the x-axis. In conclusion, consider that in the simulation, only half of the elements are really part of the simulations (and contribute to the exchange matrix PIvt)
