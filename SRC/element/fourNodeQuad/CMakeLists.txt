#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    FourNodeQuad.cpp
    FourNodeQuad3d.cpp
    EightNodeQuad.cpp
    NineNodeQuad.cpp
    PUBLIC
    FourNodeQuad.h
    FourNodeQuad3d.h
    EightNodeQuad.h
    NineNodeQuad.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
