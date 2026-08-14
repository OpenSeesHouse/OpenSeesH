#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element OBJECT)

target_sources(OPS_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      #Information.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      #Information.h
)

target_sources(OPS_Utilities
    PRIVATE
      Information.cpp
    PUBLIC
      Information.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(fourNodeQuad)
add_subdirectory(twoNodeLink)
#add_subdirectory(tetrahedron)  # pruned (Win64)
#add_subdirectory(triangle)  # pruned (Win64)
add_subdirectory(shell)
add_subdirectory(brick)


add_subdirectory(zeroLength)

add_subdirectory(truss)
#add_subdirectory(beam2d)  # pruned (Win64)
#add_subdirectory(beam3d)  # pruned (Win64)
#add_subdirectory(beamWithHinges)

add_subdirectory(dispBeamColumn)
#add_subdirectory(dispBeamColumnInt)  # pruned (Win64)
add_subdirectory(forceBeamColumn)
#add_subdirectory(mixedBeamColumn)  # pruned (Win64)
add_subdirectory(elasticBeamColumn)
#add_subdirectory(gradientInelasticBeamColumn)  # pruned (Win64)
#add_subdirectory(updatedLagrangianBeamColumn)  # pruned (Win64)
add_subdirectory(nonlinearBeamColumn)

add_subdirectory(catenaryCable)
#add_subdirectory(LHMYS)  # pruned (Win64)
#add_subdirectory(PML)  # pruned (Win64)
#add_subdirectory(RockingBC)  # pruned (Win64)
#add_subdirectory(absorbentBoundaries)  # pruned (Win64)
#add_subdirectory(adapter)  # pruned (Win64)
#add_subdirectory(componentElement)  # pruned (Win64)

add_subdirectory(frictionBearing)
add_subdirectory(elastomericBearing)
#add_subdirectory(generic)  # pruned (Win64)
add_subdirectory(joint)
#add_subdirectory(surfaceLoad)  # pruned (Win64)

#add_subdirectory(feap)  # pruned (Win64)

#add_subdirectory(PFEMElement)  # pruned (Win64)
#add_subdirectory(pipe)  # pruned (Win64)

add_subdirectory(mvlem)
#add_subdirectory(pyMacro)  # pruned (Win64)
add_subdirectory(UP-ucsd)
add_subdirectory(UWelements)
#add_subdirectory(HUelements)  # pruned (Win64)
#add_subdirectory(XMUelements)  # pruned (Win64)
#add_subdirectory(masonry)  # pruned (Win64)
#add_subdirectory(CEqElement)  # pruned (Win64)

#add_subdirectory(IGA)  # pruned (Win64)

#add_subdirectory(mefi)  # pruned (Win64)
