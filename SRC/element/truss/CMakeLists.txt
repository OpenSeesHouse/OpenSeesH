#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_truss OBJECT)

#target_sources(OPS_Element_truss
#  PRIVATE

target_sources(OPS_Element
  PRIVATE
    Truss.cpp
    TrussSection.cpp
    CorotTruss.cpp
    CorotTrussSection.cpp
  PUBLIC
    Truss.h
    TrussSection.h
    CorotTruss.h
    CorotTrussSection.h
)


#add_library(OPS_Element_InertiaTruss OBJECT)

#target_sources(OPS_Element_InertiaTruss
target_sources(OPS_Element
  PRIVATE
  PUBLIC
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

