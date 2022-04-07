module MiniTensor

using StaticArrays

const Scalar = Float64
const Vector = SVector{3, Scalar}
const Tensor = SMatrix{3, 3, Scalar}

export Scalar, Vector, Tensor

end