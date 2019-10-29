# Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems

#### *Politecnico di Milano* (ITALY)

**Author** : Annagiulia Tiozzo and Federica Laurino

**Mailto** : <annagiulia92t@gmail.com>; <federica.laurino@polimi.it>

**Date**   : July 2017


#### *Previous projects* 
**Author** : Domenico Notaro

**Mailto** : <domenico.not@gmail.com>

**Date** : March 2016

**Github Page** : https://github.com/domeniconotaro/PACS

**Author** : Stefano Brambilla

**Mailto** : <s.brambilla93@gmail.com>

**Date** : September 2016

**Github Page** : https://github.com/stefano-brambilla-853558/MANworks


-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE
- `doc/`     : Code documentation

- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources

- `config.mk`: Specify the variable GETFEM_PREFIX for GetFEM++ library

- `Doxyfile` : Instruction to build the code documentation

- `Makefile` : Instruction to install the project (see INSTALLATION)

## INSTALLATION
### Prerequisites

You need the open source finite element library "GetFEM++"

See <http://download.gna.org/getfem/html/homepage>

Version >= 4.2 is preferible

You must modify the path to the GetFEM library in `config.mk`:
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
``` 

Alternatively, at MOX cluster use the `module.sh` file:
``` 
$ source module.sh
``` 

BEWARE: 
Recall to add the library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib
```

======================

### Installation
Build the whole project with:
``` 
$ make
``` 
It first build the (static) library "libproblem3d1d" by calling
the Makefile in `include/`:
``` 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.


BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By defaul DEBUG=no.

The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include

CXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L$(GETFEM_PREFIX)/lib

LIBRARIES=-lgetfem
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

## RUN EXAMPLES
To run a specif example, go to the related subdirectory
``` 
$ cd src
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
  on the flux for the fluid problem, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel = 'DIR  DIR  DIR  DIR  DIR  DIR'

  % BC values

  BCvalue = '0.0  0.0  0.0  0.0  0.0  0.0'
  
  You can choose an arbitrary combination of
  Dirichlet-type conditions on pt and/or Robin-type conditions
  on the flux for the transport problem, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel_transp = 'MIX  MIX  MIX  MIX  MIX  MIX'

  % BC values

  BCvalue_transp = '0.0  0.0  0.0  0.0  0.0  0.0'
  
  % Coefficient for MIX condition in tissue
  
  BETAtissue_transp= 1.0E-5
  
BEWARE: All paths in file param must be ABSOLUTE

##  DEV ENVIRONMENT

GetFEM lib : 5.1
