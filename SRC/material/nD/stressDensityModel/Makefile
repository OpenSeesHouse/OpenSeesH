include ../../../../Makefile.def

# Win64 links sdmuc.lib (Fortran SDM-UC.f). Linux must build/archive the same TU.
OBJS       = stressDensity.o SDM-UC.o


all: $(OBJS)

# Optional newer 2D/3D Fortran models (not in Win64 OpenSeesPy ClCompile / sdmuc.lib)
sdm3d.o: SDM3D.f90 module_all_interfaces_3d.mod
	$(FC) $(FFLAGS) -c SDM3D.f90 -o $@

module_all_interfaces_3d.mod: MODULE_ALL_INTERFACES_3D.f90
	$(FC) $(FFLAGS) -c MODULE_ALL_INTERFACES_3D.f90

sdm2d.o: SDM2D.f90 module_all_interfaces_2d.mod
	$(FC) $(FFLAGS) -c SDM2D.f90 -o $@

module_all_interfaces_2d.mod: MODULE_ALL_INTERFACES_2D.f90
	$(FC) $(FFLAGS) -c MODULE_ALL_INTERFACES_2D.f90

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.mod

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
