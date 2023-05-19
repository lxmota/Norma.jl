@testset "schwarz-contact-dynamic-bars-hex8" begin
    cp("../examples/separate/bars/bars.yaml", "bars.yaml", force=true)
    cp("../examples/separate/bars/bar-1.yaml", "bar-1.yaml", force=true)
    cp("../examples/separate/bars/bar-2.yaml", "bar-2.yaml", force=true)
    cp("../examples/separate/bars/bar-1.g", "bar-1.g", force=true)
    cp("../examples/separate/bars/bar-2.g", "bar-2.g", force=true)
    sim = Norma.run("bars.yaml")
    subsims = sim.subsims
    model_fine = subsims[1].model
    model_coarse = subsims[2].model
    rm("bars.yaml")
    rm("bar-1.yaml")
    rm("bar-2.yaml")
    rm("bar-1.g")
    rm("bar-2.g")
    rm("bar-1.e")
    rm("bar-2.e")
    dt = sim.schwarz_controller.time_step
    t0 = sim.schwarz_controller.initial_time
    tend = sim.schwarz_controller.final_time
    nt = round((tend - t0) / dt + 1)
    nt = Int(nt)
    time = LinRange(t0, tend, nt)

    #  Analytical solution
    g = 1.e-4
    A = 1.0e-06
    rho = 1000
    E = 1.0e+09
    L1 = 2.e-3
    M1 = rho * A * L1
    V1 = 100
    t_imp = t0 + g / V1
    t_rel = t_imp + 2 * L1 * sqrt(rho / E)
    nt_fine = Int(1e+04 + 1)
    time_fine = LinRange(t0, tend, nt_fine)
    c1 = findall(time_fine .>= t_imp)
    c2 = findall(time_fine .<= t_rel)
    c = intersect(c1, c2)
    b  = findall(time_fine .< t_imp)
    a  = findall(time_fine .> t_rel) 
    contact = zeros(Bool, nt_fine)
    contact[c] .= true
    contact_force = zeros(nt_fine, 1)
    contact_force[c] .= V1 * sqrt(E * rho) * A
    velocity = zeros(nt_fine, 1)
    velocity[b] .= V1
    velocity[c] .= 0.
    velocity[a] .= -V1
    position = zeros(nt_fine, 1)
    position[b] .= V1 .* (time_fine[b] .- t0) .- g
    position[c] .= 0.
    position[a] .= -V1 .* (time_fine[a] .- t_rel)
    kinetic_energy = zeros(nt_fine, 1)
    t1 = t_imp + L1 * sqrt(rho / E)
    t2 = findall(time_fine .<= t1)
    t3 = findall(time_fine .>= t1)
    t = intersect(c1, t2)
    tt = intersect(c2, t3)
    ke_start = 0.5 * rho * A * L1 * V1 * V1
    kinetic_energy[b] .= ke_start
    kinetic_energy[t] .= ke_start .- 0.5 .* sqrt(rho * E) .* A .* V1 .* V1 .* (time_fine[t] .- t_imp)
    kinetic_energy[tt] .= 0.5 .* sqrt(rho * E) .* A .* V1 .* V1 .* (time_fine[tt] .- t1)
    kinetic_energy[a] .= ke_start
    potential_energy = zeros(nt_fine, 1)
    total_energy = zeros(nt_fine, 1)
    total_energy .= ke_start
    potential_energy = total_energy - kinetic_energy
    
    #plot(time_fine, contact, label="Analytical solution", linewidth=2)
    #plot!(time, sim.schwarz_controller.contact_hist, linecolor=:crimson, label="Schwarz", linewidth = 2)
    #title!("Contact") plot!(legend=:best) xlabel!("Time") ylabel!("Contact") xlims!(0.,8.e-6) xticks!(0:1.e-6:8.e-6) plot!(xformatter = :scientific)
    #savefig("contact.pdf")

    #plot(time_fine, contact_force, label="Analytical solution", linewidth=2)
    #plot!(time, sim.schwarz_controller.contact_hist, linecolor=:crimson, label="Schwarz", linewidth = 2)
    #title!("Contact force") plot!(legend=:best) xlabel!("Time") ylabel!("Contact force") xlims!(0.,8.e-6) xticks!(0:1.e-6:8.e-6) plot!(xformatter = :scientific)
    #savefig("contact_force.pdf")

    #plot(tV, velocity, label="Analytical solution", linewidth=2)
    #plot!(time, sim.schwarz_controller.contact_hist, linecolor=:crimson, label="Schwarz", linewidth = 2)
    #title!("Velocity") plot!(legend=:best) xlabel!("Time") ylabel!("Velocity") xlims!(0.,8.e-6) xticks!(0:1.e-6:8.e-6) plot!(xformatter = :scientific)
    #savefig("velocity.pdf")

    #plot(tV, position, label="Analytical solution", linewidth=2)
    #plot!(time, sim.schwarz_controller.contact_hist, linecolor=:crimson, label="Schwarz", linewidth = 2)
    #title!("Position") plot!(legend=:best) xlabel!("Time") ylabel!("Position") xlims!(0.,8.e-6) xticks!(0:1.e-6:8.e-6) plot!(xformatter = :scientific)
    #savefig("position.pdf")
    
end