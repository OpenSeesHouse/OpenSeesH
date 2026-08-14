#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    DoubleMembranePlateFiberSection.cpp
    ElasticMembranePlateSection.cpp
    ElasticPlateSection.cpp
    ElasticSection2d.cpp
    ElasticSection3d.cpp
    FiberSection2d.cpp
    FiberSection2dThermal.cpp
    FiberSection3d.cpp
    FiberSection3dThermal.cpp
    FiberSectionGJ.cpp
    FiberSectionGJThermal.cpp    
    LayeredShellFiberSection.cpp
    LayeredShellFiberSectionThermal.cpp
    MembranePlateFiberSection.cpp
    MembranePlateFiberSectionThermal.cpp
    NDFiberSection2d.cpp
    NDFiberSection3d.cpp
    SectionAggregator.cpp
    SectionForceDeformation.cpp
    TimoshenkoSection3d.cpp
    PUBLIC
    DoubleMembranePlateFiberSection.h
    ElasticMembranePlateSection.h
    ElasticPlateSection.h
    ElasticSection2d.h 
    ElasticSection3d.h
    FiberSection2d.h
    FiberSection2dThermal.h
    FiberSection3d.h
    FiberSection3dThermal.h
    FiberSectionGJ.h
    FiberSectionGJThermal.h
    LayeredShellFiberSection.h
    LayeredShellFiberSectionThermal.h
    MembranePlateFiberSection.h
    MembranePlateFiberSectionThermal.h
    NDFiberSection2d.h
    NDFiberSection3d.h
    SectionAggregator.h
    SectionForceDeformation.h
    TimoshenkoSection3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(fiber)
add_subdirectory(repres)
add_subdirectory(yieldSurface)
add_subdirectory(integration)
add_subdirectory(LayeredMembraneSection)
