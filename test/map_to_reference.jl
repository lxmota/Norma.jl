using LinearAlgebra
vertices= [-1  1  1 -1 -1  1  1 -1;
           -1 -1  1  1 -1 -1  1  1;
           -1 -1 -1 -1  1  1  1  1] * 0.5
element_type = "HEX8"
x1 = zeros(3)
ξ1 = Norma.map_to_reference(element_type, vertices, x1)
x2 = 0.5 * ones(3)
ξ2 = Norma.map_to_reference(element_type, vertices, x2)
x3 = -0.5 * ones(3)
ξ3 = Norma.map_to_reference(element_type, vertices, x3)
@testset "map_to_reference" begin
    @test norm(ξ1) ≈ 0.0 atol = 1.0e-08
    @test norm(ξ2 - ones(3)) ≈ 0.0 atol = 1.0e-08
    @test norm(ξ3 + ones(3)) ≈ 0.0 atol = 1.0e-08
end