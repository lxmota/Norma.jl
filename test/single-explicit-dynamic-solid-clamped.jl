cp("../examples/single/explicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml", force=true)
cp("../examples/single/explicit-dynamic-solid/clamped/clamped.g", "clamped.g", force=true)
integrator, solver, model = Norma.run("clamped.yaml")
rm("clamped.yaml")
rm("clamped.g")
rm("clamped.e")
max_disp = maximum_components(integrator.displacement)
max_velo = maximum_components(integrator.velocity)
max_acce = maximum_components(integrator.acceleration)
min_disp = minimum_components(integrator.displacement)
min_velo = minimum_components(integrator.velocity)
min_acce = minimum_components(integrator.acceleration)
@testset "single-explicit-dynamic-solid-clamped" begin
    @test max_disp[1] ≈ 0.0 atol = 1.0e-06
    @test max_disp[2] ≈ 0.0 atol = 1.0e-06
    @test max_disp[3] ≈ 0.00882497 rtol = 1.0e-04
    @test max_velo[1] ≈ 0.0 atol = 1.0e-06
    @test max_velo[2] ≈ 0.0 atol = 1.0e-06
    @test max_velo[3] ≈ 98.7857 rtol = 1.0e-03
    @test max_acce[1] ≈ 0.0 atol = 1.0e-06
    @test max_acce[2] ≈ 0.0 atol = 1.0e-06
    @test max_acce[3] ≈ 7.95633e6 rtol = 1.0e-03
    @test min_disp[1] ≈ 0.0 atol = 1.0e-06
    @test min_disp[2] ≈ 0.0 atol = 1.0e-06
    @test min_disp[3] ≈ 0.0 atol = 1.0e-06
    @test min_velo[1] ≈ 0.0 atol = 1.0e-06
    @test min_velo[2] ≈ 0.0 atol = 1.0e-06
    @test min_velo[3] ≈ -220.624 rtol = 1.0e-03
    @test min_acce[1] ≈ 0.0 atol = 1.0e-06
    @test min_acce[2] ≈ 0.0 atol = 1.0e-06
    @test min_acce[3] ≈ -1.65468e7 rtol = 1.0e-05
end
