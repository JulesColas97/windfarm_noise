# make : compilatioi
# make omp : compilation avec -fopenmp
# possibilite de decommenter lignes pour ajouter d'autres options de compilation a FCFLAGS

 SHELL = /bin/bash

chemin_modules = './'

FC = gfortran
IFC = ifort

# Warnings flags
#FCFLAGS = -Waliasing
#FCFLAGS += -Wampersand
#FCFLAGS += -Wconversion
#FCFLAGS += -Wsurprising
#FCFLAGS += -Wc-binding-type
#FCFLAGS += -Wintrinsics-std
#FCFLAGS += -Wtabs -Wall
#FCFLAGS += -Wintrinsic-shadow
#FCFLAGS += -Wline-truncation
#FCFLAGS += -Wtarget-lifetime
#FCFLAGS += -Wreal-q-constant
#FCFLAGS += std=f2003
#FCFLAGS += -fimplicit-none
FCFLAGS += -fcheck=all -fbacktrace
FCFLAGS += -O3 -fdefault-real-8  #-finit-real=nan

IFCFLAGS += -traceback -g -r8 -O3 -ipo  -fpe0  -openmp -mcmodel=medium -shared-intel

EXE = interpo
SRC = main.f90

#all: compil simu

omp : FCFLAGS += -fopenmp
omp : compil

compil :
	#Compilation des modules, l'ordre est important
	$(FC) -c $(FCFLAGS) $(chemin_modules)/bspline_kinds_module.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/bspline_sub_module.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/interpo.f90

	$(FC) -c $(FCFLAGS) $(SRC)
	
	#Compilation du programme
	# -I donne un repertoire ou chercher les .mod
#	$(FC) $(SRC)  -I$(chemin_dossier_src_comp) $(chemin_dossier_src_comp)/*.o -o  $(EXE) $(FCFLAGS)
	$(FC) $(FCFLAGS) ./*.o -o $(EXE)

clean :
	rm interpo
	rm *.o
	rm *.mod