#path to getfem library
GETFEM_PREFIX=$(mkGetfemInc)/../
#path to problem3d1d
PROBLEM_FLUID_PREFIX= /home/federica/Documenti/Software/fede_MANworks/fluid_ht_curvature/
#flag (yes/no) for debug: DEBUG=no for optimized installation
DEBUG= no
#flag (yes/no) for verbose: VERBOSE=yes prints comments at runtime
VERBOSE = no

ifeq ($(WITH_SAMG),1)
CXXFLAGS += -I${SAMG}/ -DWITH_SAMG
CXXFLAGS += -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS
LDFLAGS += -L/opt/lib/samg/ -Wl,-rpath=/opt/lib/samg/
LIBRARIES += -lamg -liomp5
endif

# getfem
CXXFLAGS+=$(shell getfem-config --cflags)
LDFLAGS+=$(shell getfem-config --libs)  

# superlu
# CXXFLAGS+=-DGMM_USES_SUPERLU -I$(mkSuperluInc)
# LDFLAGS+=-L$(mkSuperluLib) -lsuperlu
LDFLAGS+=-L$(mkQhullLib)

CXXFLAGS+=-std=c++14
CXX=g++
