#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# TODO: Create new dispBeamColumn target
target_sources(OPS_Element
  PRIVATE
    DispBeamColumn2d.cpp
    DispBeamColumn3d.cpp
    DispBeamColumn2dThermal.cpp
    DispBeamColumn3dThermal.cpp
    TimoshenkoBeamColumn2d.cpp
    TimoshenkoBeamColumn3d.cpp
  PUBLIC
    DispBeamColumn2d.h
    DispBeamColumn3d.h
    DispBeamColumn2dThermal.h
    DispBeamColumn3dThermal.h
    TimoshenkoBeamColumn2d.h
    TimoshenkoBeamColumn3d.h    
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
