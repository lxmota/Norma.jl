module MiniTensor
export MTScalar, MTVector, MTTensor

using StaticArrays

const MTScalar = Float64
const MTVector = SVector{3, MTScalar}
const MTTensor = SMatrix{3, 3, MTScalar}

end