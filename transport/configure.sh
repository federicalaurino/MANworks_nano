#Run this configure from a MOX module-based computer. 
source /u/sw/etc/profile
module load gcc-glibc/5
module load getfem
module load qhull
module load boost 

#Otherwise, set manually the paths to GetFEM, Boost and Qhull libraries:

#export mkGetfemInc=/vagrant/lib/getfem/include/
#export mkGetfemLib=/vagrant/lib/getfem/lib/
#export mkBoostInc=/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.63.0/include/
#export mkBoostLib=/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/boost/1.63.0/lib/
#export mkQhullLib=/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/qhull/2015.2/lib/
#export PATH=${mkGetfemInc}/../bin/:$PATH
#export LD_LIBRARY_PATH=${mkQhullLib}:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=${mkGetfemLib}:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=${mkBoostLib}:${LD_LIBRARY_PATH}

export mkVtkInc=/opt/VTK/include
export mkVtkLib=/opt/VTK/lib
export mkVtkHome=/opt/VTK

# SAMG libraries installed: optimized solver
export SAMG=/opt/lib/samg
export LD_LIBRARY_PATH=$SAMG:$LD_LIBRARY_PATH
export SVD_LICENSE_FILE=@nisserver.mate.polimi.iexport SVD_LICENSE_FILE=@nisserver.mate.polimi.it
# maximum number of threads
export OMP_NUM_THREADS=1
export WITH_SAMG=1
