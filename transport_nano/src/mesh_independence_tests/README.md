######################################################################################

# 	A HIERARCHIHAL MULTISCALE MODEL FOR PREDICTING
# 	THE VASCULAR BEHAVIOR OF BLOOD-BORNE NANOMEDICINES
# 	F. Laurino, A. Coclite, A. Tiozzo, P. Decuzzi, P. Zunino

######################################################################################


Here we present the mesh independent test for the Supplementaty
Material of the paper (Fig S1)

In vtk_mesh_fine: 

nb of subdivisions tissue= [44, 20, 24]
total nb of thetraedra = 190080

network mesh files = ../meshes/ratt93b
nb of subdivisions per branch = 20 
total nb of segments = 560


In vtk_mesh_fine_fine: 

nb of subdivisions tissue= [60, 40, 30]
total nb of thetraedra = 432000

network mesh files = ../meshes/ratt93b_fine
nb of subdivisions per branch = 40
total nb of segments = 1120
