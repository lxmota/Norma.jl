#
#
module LCM
export Neohookean, SaintVenant_Kirchhoff
export constitutive, isoparametric
include("operators.jl")
include("constitutive.jl")
include("interpolation.jl")
end
