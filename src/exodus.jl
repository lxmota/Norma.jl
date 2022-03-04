#
# Prepare the Exodus' Python interface for use in Julia.
#

# Point to a working Python installation. Normally this is not needed, but on
# MacOS arm64 the default Conda version seems to be broken.
ENV["PYTHON"]="/Applications/Xcode.app/Contents/Developer/usr/bin/python3"
using PyCall

# Point to a working SEACAS installation that contains exodus.py, the Exodus'
# Python interface.
push!(pyimport("sys")."path", "/Users/amota/archie/install/seacas/lib")
Exodus = pyimport("exodus")
