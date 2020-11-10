# Makefile
# 
# Usage: Initialize the houdini environment by running "source houdini_setup" in your Houdini
# directory (e.g. /opt/hfs16.5) and then run: make install
#
# Note: Eigen3 is required. Download and install it here:
#       https://eigen.tuxfamily.org/dox/index.html
#

#name of the .so that will be generated and loaded into Houdini
#this file must be inside your houdini dso folder (e.g. ~/houdini16.5/dso) to run
DSONAME = SOP_Flipsim.so 

#source files to make
SOURCES = \
    ./SOP_Flipsim.C \
    ./macgrid.C \
    ./simulator.C

#change this to your houdini dso folder
INSTDIR = ~/houdini18.0

# Additional include directories.
#INCDIRS = \
#    -I

# Additional library directories.
#LIBDIRS = \
#    -L

#set optimization flags here
#OPTIMIZER = -g

###############################################################################
# For GNU make, use this line:
#      include $(HFS)/toolkit/makefiles/Makefile.gnu
# For Microsoft Visual Studio's nmake use this line instead
#      !INCLUDE $(HFS)/toolkit/makefiles/Makefile.nmake
#
include $(HFS)/toolkit/makefiles/Makefile.gnu
