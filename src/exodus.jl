#
# Prepare the Exodus' Python interface for use in Julia.
#

# Point to a working Python installation. Normally this is not needed, but on
# MacOS arm64 the default Conda version seems to be broken.
# ENV["PYTHON"]="/Applications/Xcode.app/Contents/Developer/usr/bin/python3"
using PyCall

function exodus_module()
    # Point to a working SEACAS installation that contains exodus.py, the Exodus'
    # Python interface.
    homepath = get(ENV, "HOME", "")
    exoduspath = homepath * "/exodus/archie/install/seacas/lib"
    push!(pyimport("sys")."path", exoduspath)
    pyimport("exodus")
end
