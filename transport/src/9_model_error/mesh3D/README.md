We provide in this folder the .geo files that can build the conforming mesh in msh format.

Use conforming_mesh.geo for build a generic cubic domain that contains a cylinder with the axis parallel to the x axis. You can change the dimension of the cube, the radius and the center of the cylinder. Also, you can choose the discretation step. 
Use conforming_mesh_adaptR.geo if also you want to refine the elements around the cylinder.

All the other files are msh file containing a (-0.51, 0.51)x(-1,1)^2 domain, with a cylinder center in (x, 0,0). If non specified, the radius is 0.25. The number after conforming_mesh (e.g. conforming_mesh8) is Num_el. The 051 remembers that the x-domain is prolungated til 0.51 (in order to have the cylinder completely embedded). 
NB: "Num_el" is the number of elements you want to be in every unit of length. That is, if you set Num_el = 16, in our y-z plane (that is (-1,1)^2) along the borders you will have exactly 32 elements. 
NB: "adaptation" is the ratio between the elements around the cilynder anhd the elements on the boundary of Omega. Therefore, if adaption=0.1, around the cylinder the elements will be 10 times smaller than on the border.
