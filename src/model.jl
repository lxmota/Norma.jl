include("field.jl")

abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

struct StaticSolid <: SolidMechanics
    reference::VectorField
    current::VectorField
end

struct DynamicSolid <: SolidMechanics
    reference::VectorField
    current::VectorField
    velocity::VectorField
    acceleration::VectorField
end

struct StaticHeat <: HeatConduction
    reference::VectorField
    temperature::ScalarField
end

struct DynamicHeat <: HeatConduction
    reference::VectorField
    temperature::ScalarField
    rate::ScalarField
end