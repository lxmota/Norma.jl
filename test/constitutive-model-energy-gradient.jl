using Random

include("../src/minitensor.jl")
include("../src/constitutive_def.jl")
include("../src/constitutive.jl")

function finite_difference(material::Solid, F::Matrix{Float64}, dF::Matrix{Float64}, h::Float64)
    return (constitutive(material, F - 2*h*dF)[1] - 8*constitutive(material, F - h*dF)[1] + 8*constitutive(material, F + h*dF)[1] - constitutive(material, F + 2*h*dF)[1]) / (12*h)
end


@testset "constitutive-model-energy-gradient" begin
    Random.seed!(0)

    base_params = Dict{Any, Any}(
          "elastic modulus" => 1.0,
          "Poisson's ratio" => 0.3,
          "density" => 1.0,
    )
    sh_params = merge(base_params, Dict{Any, Any}(
          "m" => 2,
          "n" => 2,
         ))
    models = [
              Neohookean(base_params),
              SaintVenant_Kirchhoff(base_params),
              SethHill(sh_params),
             ]
    F_n = 10
    dF_n = 10

    for model in models
        for _ in 1:F_n
            F = Random.randn(3, 3) * 0.1 + I
            dWdF = constitutive(model, F)[2]
            for _ in 1:dF_n
                dF = Random.randn(3, 3) * 0.1
                dWdF_fd = finite_difference(model, F, dF, 1e-6)
                @test isapprox(sum(dWdF .* dF), dWdF_fd, atol=1e-6)
            end
        end
    end
end