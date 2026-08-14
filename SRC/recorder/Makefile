include ../../Makefile.def

ifeq ($(GRAPHICS), NONE)
GRAPHIC_OBJECTS =
else
GRAPHIC_OBJECTS = FilePlotter.o
endif

OBJS       = Recorder.o \
	ElementRecorder.o \
	NodeRecorder.o \
	EnvelopeElementRecorder.o \
	EnvelopeNodeRecorder.o \
	DriftRecorder.o \
	EnvelopeDriftRecorder.o \
	$(GRAPHIC_OBJECTS) \
	ResidElementRecorder.o \
	ResidNodeRecorder.o \
	ResidDriftRecorder.o \
	ConditionalElementRecorder.o \
	ConditionalNodeRecorder.o \
	ConditionalDriftRecorder.o


# Compilation control

all:      $(OBJS) $(GRAPHIC_OBJECTS)
	@$(CD) $(FE)/recorder/response; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/recorder/response; $(MAKE) wipe;

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
