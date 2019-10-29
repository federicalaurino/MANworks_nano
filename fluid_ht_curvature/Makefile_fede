# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : November 2015
# ====================================================================

.PHONY: all doc clean distclean library fast

all: library
	$(MAKE) -C src/1_curved_singlebranch
	$(MAKE) -C src/2_bifurcation
	$(MAKE) -C src/3_anastomosis
	$(MAKE) -C src/4_Voronoi_Network

library: 
	$(MAKE) -C include

doc:
	install -d doc
	doxygen Doxyfile
	
pdf: doc
	$(MAKE) -C doc/latex pdf

clean:
	$(RM) -r *~ *.log
	$(MAKE) -C include clean
	$(MAKE) -C src/1_curved_singlebranch clean
	$(MAKE) -C src/2_bifurcation clean
	$(MAKE) -C src/3_anastomosis clean
	$(MAKE) -C src/4_Voronoi_Network clean

distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_curved_singlebranch distclean
	$(MAKE) -C src/2_bifurcation distclean
	$(MAKE) -C src/3_anastomosis distclean
	$(MAKE) -C src/4_Voronoi_Network distclean

fast:
	$(MAKE) -C include fastclean
	$(MAKE) -C include
	$(MAKE) -C src/1_curved_singlebranch  fastclean
	$(MAKE) -C src/2_bifurcation fastclean
	$(MAKE) -C src/3_anastomosis fastclean
	$(MAKE) -C src/4_Voronoi_Network fastclean

