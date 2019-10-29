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


# The SAMG libraries are not installed on this computer, non optimized solver.
export WITH_SAMG=0
