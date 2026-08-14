#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    Domain.cpp
    DomainModalProperties.cpp
  PUBLIC
    Domain.h
    DomainModalProperties.h
)

add_subdirectory(single)
add_subdirectory(partitioned)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

