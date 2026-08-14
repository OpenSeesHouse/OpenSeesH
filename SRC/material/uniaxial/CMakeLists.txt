#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
      AxialSp.cpp
      AxialSp.cpp
      AxialSpHD.cpp
      Bilin.cpp
      ConcreteCM.cpp
      ElasticMaterial.cpp
      ElasticMultiLinear.cpp
      ElasticPPMaterial.cpp
      IMKBilin.cpp
      IMKPeakOriented.cpp
      IMKPinching.cpp
      MaterialState.cpp
      ParallelMaterial.cpp
      UniaxialMaterial.cpp
    PUBLIC
      AxialSp.h
      AxialSp.h
      AxialSpHD.h
      Bilin.h
      ConcreteCM.h
      ElasticMaterial.h
      ElasticMultiLinear.h
      ElasticPPMaterial.h
      IMKBilin.h
      IMKPeakOriented.h
      IMKPinching.h
      MaterialState.h
      ParallelMaterial.h
      UniaxialMaterial.h
)


target_sources(OPS_Material
    PRIVATE
      BarSlipMaterial.cpp
      Bond_SP01.cpp
      CableMaterial.cpp
      Concrete01.cpp
      Concrete02.cpp
      Concrete02Thermal.cpp
      Concrete04.cpp
      Concrete05.cpp
      Concrete06.cpp
      Concrete07.cpp
      ConcreteECThermal.cpp
      DamperMaterial.cpp
      ENTMaterial.cpp
      EPPGapMaterial.cpp
      ElasticMaterialThermal.cpp
      FatigueMaterial.cpp
      #FedeasMaterial.cpp
      HardeningMaterial.cpp
      #HardeningMaterial2.cpp
      HyperbolicGapMaterial.cpp
      HystereticMaterial.cpp
      HystereticSmooth.cpp
      InitStressMaterial.cpp
      KikuchiAikenHDR.cpp
      KikuchiAikenLRB.cpp
      #Maxwell.cpp                      # pruned (Win64)
      MinMaxMaterial.cpp
      ModIMKPeakOriented.cpp
      ModIMKPeakOriented02.cpp
      ModIMKPinching.cpp
      ModIMKPinching02.cpp
      PathIndependentMaterial.cpp
      PenaltyMaterial.cpp      
      Pinching4Material.cpp
      RambergOsgoodSteel.cpp
      ReinforcingSteel.cpp
      SPSW02.cpp
      SelfCenteringMaterial.cpp
      SeriesMaterial.cpp
      ShearPanelMaterial.cpp
      StainlessECThermal.cpp
      Steel01.cpp
      Steel01Thermal.cpp
      Steel02.cpp
      Steel02Thermal.cpp
      Steel4.cpp
      SteelECThermal.cpp
      UniaxialJ2Plasticity.cpp
      ViscousDamper.cpp
      ViscousMaterial.cpp
      pyUCLA.cpp
      SMAMaterial.cpp
      Steel05.cpp
      SmoothIMK.cpp
      ConfinedConcrete.cpp
      Steel_01.cpp
       

    PUBLIC
      BarSlipMaterial.h
      Bond_SP01.h
      CableMaterial.h
      Concrete01.h
      Concrete02.h
      Concrete02Thermal.h
      Concrete04.h
      Concrete05.h
      Concrete06.h
      Concrete07.h
      ConcreteECThermal.h
      ConfinedConcrete.h
      DamperMaterial.h
      ENTMaterial.h
      EPPGapMaterial.h
      ElasticMaterialThermal.h
      FatigueMaterial.h
      #FedeasMaterial.h
      HardeningMaterial.h
      #HardeningMaterial2.h
      HyperbolicGapMaterial.h
      HystereticMaterial.h
      HystereticSmooth.h
      InitStressMaterial.h
      KikuchiAikenHDR.h
      KikuchiAikenLRB.h
      MinMaxMaterial.h
      ModIMKPeakOriented.h
      ModIMKPeakOriented02.h
      ModIMKPinching.h
      ModIMKPinching02.h
      PathIndependentMaterial.h
      PenaltyMaterial.h
      Pinching4Material.h
      RambergOsgoodSteel.h
      ReinforcingSteel.h
      SPSW02.h
      SelfCenteringMaterial.h
      SeriesMaterial.h
      ShearPanelMaterial.h
      SmoothIMK.h
      SMAMaterial.h
      StainlessECThermal.h
      Steel01.h
      Steel01Thermal.h
      Steel02.h
      Steel02Thermal.h
      Steel4.h
      Steel05.h
      Steel_01.h
      SteelECThermal.h
      UniaxialJ2Plasticity.h
      ViscousDamper.h
      ViscousMaterial.h
      pyUCLA.h

      TriMatrix.h  # move to SRC/matrix ?
)



target_sources(OPS_Material
  PRIVATE
  PUBLIC
)

target_sources(OPS_Material
  PRIVATE
    QbSandCPT.cpp
    TzSandCPT.cpp
  PUBLIC
    QbSandCPT.h
    TzSandCPT.h
)

add_library(OPS_Material_f)

target_sources(OPS_Material_f PRIVATE
      DoddRestrepo.f
      STEELDR.f
      c14-SK-M.f
)

target_link_libraries(OPS_Material PUBLIC OPS_Material_f)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#add_subdirectory(backbone)  # pruned (Win64)
#add_subdirectory(stiffness)  # pruned (Win64)
#add_subdirectory(strength)  # pruned (Win64)
#add_subdirectory(unloading)  # pruned (Win64)
#add_subdirectory(limitState)  # pruned (Win64)

add_subdirectory(PY)
#add_subdirectory(snap)  # pruned (Win64)
#add_subdirectory(fedeas)  # pruned (Win64)
#add_subdirectory(drain)  # pruned (Win64)



