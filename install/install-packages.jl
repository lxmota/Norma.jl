using Pkg

# Point to a working Python installation. Normally this is not needed, but on
# MacOS arm64 the default Conda version seems to be broken.
# ENV["PYTHON"]="/Applications/Xcode.app/Contents/Developer/usr/bin/python3"

Pkg.add("DelimitedFiles")
Pkg.add("Formatting")
Pkg.add("LinearAlgebra")
Pkg.add("PyCall")
Pkg.add("SparseArrays")
Pkg.add("StaticArrays")
Pkg.add("Symbolics")