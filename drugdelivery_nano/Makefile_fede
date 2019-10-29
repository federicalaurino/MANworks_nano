# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright D. Brambilla 2016
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com>
#   DATE        : Settembre 2016
# ====================================================================

.PHONY: all doc clean distclean library

all: library
	$(MAKE) -C src/1_singlebranch_transp 
	$(MAKE) -C src/2_Ybifurcation_transp 
	$(MAKE) -C src/3_network_transp 
	$(MAKE) -C src/4_curved_singlebranch 
	$(MAKE) -C src/5_curved_bifurcation
	$(MAKE) -C src/6_anastomosis 
	$(MAKE) -C src/7_Voronoi_Network
	$(MAKE) -C src/8_convergence_error
	$(MAKE) -C src/9_model_error

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
	$(MAKE) -C src/1_singlebranch_transp clean
	$(MAKE) -C src/2_Ybifurcation_transp clean
	$(MAKE) -C src/3_network_transp clean
	$(MAKE) -C src/4_curved_singlebranch clean
	$(MAKE) -C src/5_curved_bifurcation clean
	$(MAKE) -C src/6_anastomosis  clean
	$(MAKE) -C src/7_Voronoi_Network clean
	$(MAKE) -C src/8_convergence_error clean
	$(MAKE) -C src/9_model_error clean
	
distclean: clean
	$(RM) -r doc/*
	$(MAKE) -C include distclean
	$(MAKE) -C src/1_singlebranch_transp distclean
	$(MAKE) -C src/2_Ybifurcation_transp distclean
	$(MAKE) -C src/3_network_transp distclean
	$(MAKE) -C src/4_curved_singlebranch distclean
	$(MAKE) -C src/5_curved_bifurcation distclean
	$(MAKE) -C src/6_anastomosis  distclean
	$(MAKE) -C src/7_Voronoi_Network distclean
	$(MAKE) -C src/8_convergence_error distclean
	$(MAKE) -C src/9_model_error distclean
