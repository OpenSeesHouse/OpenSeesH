include ../../../Makefile.def

OBJS       = UniaxialMaterial.o \
	ElasticMaterial.o \
	Bilin.o \
	IMKBilin.o \
	IMKPeakOriented.o \
	IMKPinching.o \
	ConcreteCM.o \
	ElasticMultiLinear.o \
	ElasticPPMaterial.o \
	ParallelMaterial.o \
	AxialSp.o \
	AxialSpHD.o \
	ViscousDamper.o \
	HardeningMaterial.o \
	Concrete01.o \
	Concrete02Thermal.o \
	Steel02Thermal.o \
	Steel01Thermal.o \
	Concrete04.o \
	Concrete06.o \
	Concrete07.o \
	RambergOsgoodSteel.o \
	Steel01.o \
	HystereticMaterial.o \
	EPPGapMaterial.o \
	ViscousMaterial.o \
	SeriesMaterial.o \
	InitStressMaterial.o \
	MinMaxMaterial.o \
	PathIndependentMaterial.o \
	CableMaterial.o \
	ENTMaterial.o \
	Pinching4Material.o \
	BarSlipMaterial.o \
	FatigueMaterial.o \
	Steel02.o \
	Concrete02.o \
	Bond_SP01.o \
	UniaxialJ2Plasticity.o \
	HyperbolicGapMaterial.o \
	ReinforcingSteel.o \
	PenaltyMaterial.o \
	KikuchiAikenHDR.o \
	KikuchiAikenLRB.o \
	Concrete05.o \
	ShearPanelMaterial.o \
	SelfCenteringMaterial.o \
	ModIMKPeakOriented.o \
	ModIMKPeakOriented02.o \
	ModIMKPinching.o \
	ModIMKPinching02.o \
	pyUCLA.o \
	Steel4.o \
	ConcreteECThermal.o \
	ElasticMaterialThermal.o \
	StainlessECThermal.o \
	SteelECThermal.o \
	DamperMaterial.o \
	MaterialState.o \
	SPSW02.o \
	SMAMaterial.o \
	Steel05.o \
	HystereticSmooth.o \
	QbSandCPT.o \
	TzSandCPT.o \
	SmoothIMK.o \
	ConfinedConcrete.o \
	Steel_01.o


all:         $(OBJS)
	@$(CD) $(FE)/material/uniaxial/PY; $(MAKE);
#	@$(CD) $(FE)/material/uniaxial/fedeas; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/drain; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/snap; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/backbone; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/stiffness; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/strength; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/unloading; $(MAKE);  # pruned (Win64)
#	@$(CD) $(FE)/material/uniaxial/limitState; $(MAKE);  # pruned (Win64)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/uniaxial/PY; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/fedeas; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/drain; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/snap; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/backbone; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/stiffness; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/strength; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/unloading; $(MAKE) wipe;
#	@$(CD) $(FE)/material/uniaxial/limitState; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
