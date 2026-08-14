#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        Joint2D.cpp
        Joint3D.cpp
        MP_Joint2D.cpp
        MP_Joint3D.cpp
        #TclBeamColumnJointCommand.cpp
        #TclJoint2dCommand.cpp
        #TclJoint3dCommand.cpp
    PUBLIC
        Joint2D.h
        Joint3D.h
        MP_Joint2D.h
        MP_Joint3D.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

