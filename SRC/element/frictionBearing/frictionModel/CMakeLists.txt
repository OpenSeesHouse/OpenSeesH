#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        Coulomb.cpp
        FrictionModel.cpp
        FrictionResponse.cpp
    PUBLIC
        Coulomb.h
        FrictionModel.h
        FrictionResponse.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

