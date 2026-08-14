#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
    ZeroLength.cpp
    ZeroLengthND.cpp
    ZeroLengthSection.cpp
    PUBLIC
    ZeroLength.h
    ZeroLengthND.h
    ZeroLengthSection.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
