include ../../Makefile.def

OBJS       = Element.o Information.o ElementalLoad.o

# Compilation control — recurse only into Win64-allowlisted element trees

all:         $(OBJS)
	@$(CD) $(FE)/element/UP-ucsd; $(MAKE);
	@$(CD) $(FE)/element/UWelements; $(MAKE);
	@$(CD) $(FE)/element/brick; $(MAKE);
	@$(CD) $(FE)/element/catenaryCable; $(MAKE);
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/elastomericBearing; $(MAKE);
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/fourNodeQuad; $(MAKE);
	@$(CD) $(FE)/element/frictionBearing; $(MAKE);
	@$(CD) $(FE)/element/joint; $(MAKE);
	@$(CD) $(FE)/element/mvlem; $(MAKE);
	@$(CD) $(FE)/element/nonlinearBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/shell; $(MAKE);
	@$(CD) $(FE)/element/surfaceLoad; $(MAKE);
	@$(CD) $(FE)/element/truss; $(MAKE);
	@$(CD) $(FE)/element/twoNodeLink; $(MAKE);
	@$(CD) $(FE)/element/zeroLength; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) Makefile.bak $(OBJS) *.o *~ #*# core

spotless: clean
	@$(CD) $(FE)/element/UP-ucsd; $(MAKE) wipe;
	@$(CD) $(FE)/element/UWelements; $(MAKE) wipe;
	@$(CD) $(FE)/element/brick; $(MAKE) wipe;
	@$(CD) $(FE)/element/catenaryCable; $(MAKE) wipe;
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/elastomericBearing; $(MAKE) wipe;
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/fourNodeQuad; $(MAKE) wipe;
	@$(CD) $(FE)/element/frictionBearing; $(MAKE) wipe;
	@$(CD) $(FE)/element/joint; $(MAKE) wipe;
	@$(CD) $(FE)/element/mvlem; $(MAKE) wipe;
	@$(CD) $(FE)/element/nonlinearBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/shell; $(MAKE) wipe;
	@$(CD) $(FE)/element/surfaceLoad; $(MAKE) wipe;
	@$(CD) $(FE)/element/truss; $(MAKE) wipe;
	@$(CD) $(FE)/element/twoNodeLink; $(MAKE) wipe;
	@$(CD) $(FE)/element/zeroLength; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
