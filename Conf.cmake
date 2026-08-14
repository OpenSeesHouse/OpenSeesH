#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#==============================================================================


#==============================================================================
#                           Select Default Executable
#==============================================================================
set(OPS_FINAL_TARGET "OpenSeesPy" CACHE STRING "OpenSees final target")

# When OFF, skip Tcl find/link/install and Tcl interpreter targets (Python-only).
option(OPS_WITH_TCL "Build with Tcl interpreter support" OFF)

#==============================================================================
#                            Basic Switches
#==============================================================================

option(OPS_Use_Dev_Directories
  "Include files in DEVELOPER directory"                   OFF)

option(FMK
  "Special FMK Code"                                       OFF)

set(OPS_Use_Graphics_Option
  None
  # Base
  # OpenGL
)

