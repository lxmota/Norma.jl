# This is needed for now on Mac but let's hope it goes away.
ENV["JULIA_PKG_IGNORE_HASHES"] = 1

using Pkg

Pkg.add("DelimitedFiles")
Pkg.add("Exodus")
Pkg.add("Formatting")
Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
Pkg.add("Symbolics")
Pkg.add("YAML")