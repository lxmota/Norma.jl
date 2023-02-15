cp("../examples/single/implicit-dynamic-solid/clamped/clamped.yaml", "clamped.yaml", force=true)
cp("../examples/single/implicit-dynamic-solid/clamped/clamped.g", "clamped.g", force=true)
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
@testset "single-implicit-dynamic-solid-clamped" begin
    @test max_disp[1] ≈ 0.0 atol = 1.0e-06
    @test max_disp[2] ≈ 0.0 atol = 1.0e-06
    @test max_disp[3] ≈ 0.008825600516255708 rtol = 1.0e-06
    @test max_velo[1] ≈ 0.0 atol = 1.0e-06
    @test max_velo[2] ≈ 0.0 atol = 1.0e-06
    @test max_velo[3] ≈ 98.78220185628787 rtol = 1.0e-06
    @test max_acce[1] ≈ 0.0 atol = 1.0e-06
    @test max_acce[2] ≈ 0.0 atol = 1.0e-06
    @test max_acce[3] ≈ 7.959961754746247e6 rtol = 1.0e-06
    @test min_disp[1] ≈ 0.0 atol = 1.0e-06
    @test min_disp[2] ≈ 0.0 atol = 1.0e-06
    @test min_disp[3] ≈ 0.0 atol = 1.0e-06
    @test min_velo[1] ≈ 0.0 atol = 1.0e-06
    @test min_velo[2] ≈ 0.0 atol = 1.0e-06
    @test min_velo[3] ≈ -220.6515070475273 rtol = 1.0e-06
    @test min_acce[1] ≈ 0.0 atol = 1.0e-06
    @test min_acce[2] ≈ 0.0 atol = 1.0e-06
    @test min_acce[3] ≈ -1.6561220204901198e7 rtol = 1.0e-06
end
