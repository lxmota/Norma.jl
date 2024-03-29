@testset "single-implicit-dynamic-solid-sho" begin
    cp("../examples/single/implicit-dynamic-solid/sho/sho.yaml", "sho.yaml", force=true)
    cp("../examples/single/implicit-dynamic-solid/sho/sho.g", "sho.g", force=true)
    simulation = Norma.run("sho.yaml")
    integrator = simulation.integrator
    rm("sho.yaml")
    rm("sho.g")
    rm("sho.e")
    displacement = integrator.displacement[end]
    velocity = integrator.velocity[end]
    acceleration = integrator.acceleration[end]
    @test displacement ≈ 0.0 atol = 2.0e-03
    @test velocity ≈ -1.0 rtol = 1.0e-05
    @test acceleration ≈ 0.0 atol = 2.0e-03
end
