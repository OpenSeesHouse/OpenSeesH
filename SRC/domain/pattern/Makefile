include ../../../Makefile.def

# Random-process series pruned (Win64); keep empty even if RELIABILITY=YES
RELIABILITY_OBJS =

OBJS       = EarthquakePattern.o \
	LinearSeries.o \
	LoadPattern.o \
	FireLoadPattern.o \
	LoadPatternIter.o \
	PathSeries.o \
	PathTimeSeries.o \
	PathTimeSeriesThermal.o \
	RectangularSeries.o \
	TimeSeries.o \
	TrigSeries.o \
	RampSeries.o \
	ConstantSeries.o \
	MPAccSeries.o \
	UniformExcitation.o \
	TimeSeriesIntegrator.o \
	TrapezoidalTimeSeriesIntegrator.o \
	SimpsonTimeSeriesIntegrator.o \
	MultiSupportPattern.o \
	PulseSeries.o \
	TriangleSeries.o \
	PeerMotion.o \
	PeerNGAMotion.o


# Compilation control

all:         $(OBJS)
#	@$(CD) $(FE)/domain/pattern/drm; $(MAKE);  # pruned (Win64)


# Miscellaneous
tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
#	@$(CD) $(FE)/domain/pattern/drm; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.
