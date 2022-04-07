abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

struct StaticSolid <: SolidMechanics
end

struct DynamicSolid <: SolidMechanics
end

struct StaticHeat <: HeatConduction
end

struct DynamicHeat <: HeatConduction
end