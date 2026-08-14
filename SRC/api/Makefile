include ../../Makefile.def

# packages.o: dynamic plugin loader (dlopen) — needed by ObjectBroker for both Tcl/Python
# elementAPI_TCL.o / elementAPI_Dummy.o: Tcl-only or obsolete; OPS_* live in
# SRC/interpreter/OpenSeesCommands.cpp for OpenSeesPy
OBJS       = packages.o

# Compilation control

all: 	$(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
