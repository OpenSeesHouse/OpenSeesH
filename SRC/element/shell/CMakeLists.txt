#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        R3vectors.cpp
        ShellDKGQ.cpp
        ShellDKGT.cpp
        ShellMITC4.cpp
        ShellMITC4Thermal.cpp
        ShellMITC9.cpp
        ShellNLDKGQ.cpp
        ShellNLDKGQThermal.cpp
        ShellNLDKGT.cpp
        ShellNLDKGTThermal.cpp
    PUBLIC
        R3vectors.h
        ShellDKGQ.h
        ShellDKGT.h
        ShellMITC4.h
        ShellMITC4Thermal.h
        ShellMITC9.h
        ShellNLDKGQ.h
        ShellNLDKGQThermal.h
        ShellNLDKGT.h
        ShellNLDKGTThermal.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

