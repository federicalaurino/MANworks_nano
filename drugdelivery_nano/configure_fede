#Run this configure from a MOX module-based computer. 
#source /u/sw/etc/profile
#module load gcc-glibc/5
#module load getfem
#module load qhull
#module load boost

#Otherwise, set manually the paths to GetFEM, Boost and Qhull libraries:

export mkGetfemInc=/usr/local/getfem_5_1/include/
export mkGetfemLib=/usr/local/getfem_5_1/lib
export mkBoostInc=/usr/include/boost/
export mkBoostLib=/usr/lib/x86_64-linux-gnu/
export mkQhullLib=/usr/local/lib
export PATH=${mkGetfemInc}/../bin/:$PATH
export LD_LIBRARY_PATH=${mkQhullLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkGetfemLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkBoostLib}:${LD_LIBRARY_PATH}

 
export mkVtkInc=/opt/VTK/include
export mkVtkLib=/opt/VTK/lib
export mkVtkHome=/opt/VTK

# SAMG libraries installed: optimized solver
#export SAMG=/opt/lib/samg
#export LD_LIBRARY_PATH=$SAMG:$LD_LIBRARY_PATH
#export SVD_LICENSE_FILE=@nisserver.mate.polimi.iexport SVD_LICENSE_FILE=@nisserver.mate.polimi.it
# maximum number of threads
export OMP_NUM_THREADS=1
export WITH_SAMG=0
