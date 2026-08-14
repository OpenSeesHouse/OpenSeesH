include ../../../Makefile.def

OBJS       = FourNodeQuad.o \
	FourNodeQuad3d.o \
	NineNodeQuad.o \
	EightNodeQuad.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
