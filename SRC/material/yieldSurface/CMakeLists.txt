#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Material_YieldSurface OBJECT)

target_include_directories(OPS_Material_YieldSurface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Sources pruned to match Win64 (empty object library retained for OpenSeesLIB)
#add_subdirectory(./evolution/)
#add_subdirectory(./plasticHardeningMaterial/)
#add_subdirectory(./yieldSurfaceBC/)
