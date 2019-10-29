# Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems
#### *Politecnico di Milano* (ITALY)

**Author** : Stefano Brambilla 

**Mailto** : <s.brambilla93@gmail.com>

**Date**   : March 2018

#### *Previous projects*

**Authors** :  Luca Possenti, Simone Di Gregorio, Giorgio Raimondi, Fannie Gerosa

**Mailto** : <luca.possenti@polimi.it>

**Date**   : January 2018

**Github Page** : https://github.com/lpossenti/MANworks_ht_curvature

---------------------
**Author** : Domenico Notaro 

**Mailto** : <domenico.not@gmail.com>

**Date**   : March 2016

**Github Page** : https://github.com/domeniconotaro/PACS

-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE

- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources
  - `1_singlebranch_transp/` : solve the coupling with single-vessel network
  - `2_Ybifurcation_transp/` : solve the problem with Y-shaped network
  - `3_network_transp/`  : solve te problem with a small network
  - `4_curved_singlebranch/`  : solve the problem with a curved single-vessel network
  - `5_curved_bifurcation/`  : solve the problem with curved Y-shaped network
  - `6_anastomosis/`  : solve the problem with anastomosis network
  - `7_Voronoi_network/`  : solve te problem with a Voronoi network
 
  
- `config.mk`: Specify some variables for compiler

- `configure.sh`: Load the modules from a MOX computer

- `Doxyfile` : Instruction to build the code documentation

- `Makefile` : Instruction to install the project (see INSTALLATION)

## INSTALLATION
### Prerequisites

You need the open source finite element library "GetFEM++"

See <http://download.gna.org/getfem/html/homepage>

Version >= 5.1 is necessary

BEWARE: 
Recall to add the library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib
```
======================

### Before installation

You must modify, in `config.mk`, the path to the GetFEM library and to folder containing the fluid problem:
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
PROBLEM_FLUID_PREFIX=/home/.../path/to/.../fluid_ht_curvature
``` 
In `config.mk`, you can also set the flags for optimized or debug installation, and for activating the comments at runtime:
``` 
DEBUG= no
VERBOSE= yes
``` 

Finally, you need to export to LD_LIBRARY_PATH the paths to GetFEM, Boost and Qhull libraries;
this can be done using the modules system (from a MOX computer) or setting manually the paths in the file configure.sh.
Before compiling, call:
``` 
$ source configure.sh
``` 

======================

### Installation
Build the whole project with:
``` 
$ make
``` 
It first build the (shared) library "libproblem3d1d" by calling
the Makefile in `include/`:
``` 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.

It is also possible to build a single example, e.g. "1_singlebranch_transp", with:
``` 
$ make library

$ make -C src/1_singlebranch_transp
``` 


The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include -I$(mkBoostInc) -I$(PROBLEM_FLUID_PREFIX)/include

CXXFLAGS=-std=c++11 

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L../../lib -L$(GETFEM_PREFIX)/lib -L$(PROBLEM_FLUID_PREFIX)/lib -L$(mkBoostLib) -Wl,-rpath=../../lib -Wl,-rpath=$(PROBLEM_FLUID_PREFIX)/lib 

LIBRARIES=-lgetfem -lproblem3d1d -ltransport3d1d -lutil -lboost_iostreams -lboost_system -lboost_filesystem
``` 
Recall that any macro may be overrlued by specifying it when calling 
make. Example: 
``` 
$ make CXXFLAGS+=-DSOMETHING OPTFLAGS=-g
``` 

======================

### Documentation
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince doc/latex/refman.pdf
``` 

## MAKE OPTIONS
All examples are provided with a Makefile which accepts the following
options:
-  all       : makes the example
-  clean     : as it says
-  distclean : clean and also deletes temporary file and local doc directory
Being "all" the first target of the makefile, to compile the examples is
sufficient to type make. 
In addition the external Makefile (./Makefile) has the following options:
-  doc       : produces the documentation (html, tex)
-  pdf       : produces a guide in portable format
- library    : build the library from files in include/
- static     : build a static library from the files in include/

## RUN EXAMPLES
To run a specific example, go to the related subdirectory
``` 
$ cd src/3_network_transp
``` 
Build the program
``` 
$ make
``` 
Execute the program with specific input
``` 
$ ./M3D1D input.param
``` 
Each program contains the file input.param that specifies 

- Some flags to identify the particular example
  -  TEST_PARAM = 1  # import parameters in dimensionless form
  -  VTK_EXPORT = 1  # export results in vtk format
  -  ...

- The mesh
  - For the 3D mesh you can either provide instruction to build a simple
  regular mesh (TEST_GEOMETRY = 1) or the absolute path to import a mesh
  pre-built with Gmsh (.gmsh)
  - For the 1D mesh specify the path to the file of points (.pts). All
  examples come with a possible pts file

- GetFEM++ descriptors (FEM, ...)

- Problem parameters (dimensional or dimensionless)

- Boundary conditions. You can choose an arbitrary combination of
  Dirichlet-type conditions on pt and/or Robin-type conditions
  on the flux, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel = 'DIR  DIR  DIR  DIR  DIR  DIR'

  % BC values

  BCvalue = '0.0  0.0  0.0  0.0  0.0  0.0'
  
 

##  DEV ENVIRONMENT
OS         : CentOS Linux 7 64-bit 

Processor  : Intel® Core™ i5-2310 CPU @ 2.90GHz × 4

Compiler   : g++-5.2.1

GetFEM lib : 5.2

