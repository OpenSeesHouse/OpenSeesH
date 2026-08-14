#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Recorder
    PRIVATE
      DriftRecorder.cpp
      ElementRecorder.cpp
      EnvelopeDriftRecorder.cpp
      EnvelopeElementRecorder.cpp
      EnvelopeNodeRecorder.cpp
      NodeRecorder.cpp
      Recorder.cpp
      ResidDriftRecorder.cpp
      ResidElementRecorder.cpp
      ResidNodeRecorder.cpp
      ConditionalDriftRecorder.cpp
      ConditionalElementRecorder.cpp
      ConditionalNodeRecorder.cpp
    PUBLIC
      DriftRecorder.h
      ElementRecorder.h
      EnvelopeDriftRecorder.h
      EnvelopeElementRecorder.h
      EnvelopeNodeRecorder.h
      NodeRecorder.h
      Recorder.h
      ResidDriftRecorder.h
      ResidElementRecorder.h
      ResidNodeRecorder.h
      ConditionalDriftRecorder.h
      ConditionalElementRecorder.h
      ConditionalNodeRecorder.h
)

target_sources(OPS_Graphics
    PRIVATE
      FilePlotter.cpp
    PUBLIC
      FilePlotter.h
)

add_subdirectory(response)
