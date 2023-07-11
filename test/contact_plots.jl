using CSV
using DataFrames
using YAML
using Plots
using Formatting
using LinearAlgebra

dt = sim.schwarz_controller.time_step
t0 = sim.schwarz_controller.initial_time
tend = sim.schwarz_controller.final_time
nt = round((tend - t0) / dt + 1)
nt = Int(nt)
time = LinRange(t0, tend, nt)

disp_x1 = zeros(nt)
velo_x1 = zeros(nt)
acce_x1 = zeros(nt)
potential_energy_sim = zeros(nt)
kinetic_energy_sim = zeros(nt)

for i = 1:nt
    stop = i - 1
    index_string = "-" * string(stop, pad=4)
    disp_tot = CSV.read("01-disp" * index_string *  ".csv", DataFrame,  header=false)[!,1]
    velo_tot = CSV.read("01-velo" * index_string *  ".csv", DataFrame,  header=false)[!,1]
    acce_tot = CSV.read("01-acce" * index_string *  ".csv", DataFrame,  header=false)[!,1]
    kinetic_energy_tot = CSV.read("01-kinetic" * index_string *  ".csv", DataFrame,  header=false)[!,1]
    potential_energy_tot = CSV.read("01-potential" * index_string *  ".csv", DataFrame,  header=false)[!,1]
    dof = size(disp_tot)
    disp_x1[i] = disp_tot[3*183-2]
    velo_x1[i] = velo_tot[3*183-2]
    acce_x1[i] = acce_tot[3*183-2]
    potential_energy_sim[i] = potential_energy_tot[1]
    kinetic_energy_sim[i] = kinetic_energy_tot[1]
    rm("01-disp" * index_string *  ".csv")
    rm("02-disp" * index_string *  ".csv")
    rm("01-velo" * index_string *  ".csv")
    rm("02-velo" * index_string *  ".csv")
    rm("01-acce" * index_string *  ".csv")
    rm("02-acce" * index_string *  ".csv")
    rm("01-kinetic" * index_string *  ".csv")
        m("02-potential" * index_string *  ".csv")
end

# Analytical solution
g = 1.e-4
A = 1.0e-08
rho = 1000
E = 1.0e+09
L1 =1.e-3
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
total_en = zeros(nt, 1)
total_en .= ke_start  
    
plot_figs = false

if plot_figs == true
    plot(time_fine, contact, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, sim.schwarz_controller.contact_hist, linecolor=:crimson, label="Schwarz", linewidth = 2)
    title!("Contact") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("Contact") 
    plot!(xformatter = :scientific)
    savefig("contact.pdf")

    plot(time_fine, velocity, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, velo_x1, linecolor=:crimson, label="Schwarz point 1", linewidth = 2)
    title!("Velocity") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("Velocity") 
    plot!(xformatter = :scientific)
    savefig("velocity.pdf")

    plot(time_fine, position, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, disp_x1,  linecolor=:crimson, label="Schwarz point 1" , linewidth = 2)
    title!("displacement") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("displacement")
    plot!(xformatter = :scientific)
    savefig("displacement.pdf")

    plot(time_fine, kinetic_energy, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, kinetic_energy_sim, linecolor=:crimson, label="Schwarz point 1" , linewidth = 2)
    title!("kinetic energy") 
    plot!(legend=:top) 
    xlabel!("Time") 
    ylabel!("kinetic energy") 
    plot!(xformatter = :scientific)
    savefig("kinetic.pdf")

    plot(time_fine, potential_energy, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, potential_energy_sim, linecolor=:crimson, label="Schwarz point 1" , linewidth = 2)
    title!("potential energy") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("potential energy") 
    plot!(xformatter = :scientific)
    savefig("potential.pdf")

    plot(time_fine, total_energy, linecolor=:black, label="Analytical solution", linewidth=2)
    plot!(time, (potential_energy_sim+kinetic_energy_sim), linecolor=:crimson, label="Schwarz point 1" , linewidth = 2)
    title!("total energy") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("total energy") 
    plot!(xformatter = :scientific)
    savefig("total.pdf")

    plot(time, ((potential_energy_sim+kinetic_energy_sim) - total_en) ./ total_en, linecolor=:crimson, label="Schwarz point 1" , linewidth = 2)
    title!("total energy") 
    plot!(legend=:bottom) 
    xlabel!("Time") 
    ylabel!("total energy") 
    plot!(xformatter = :scientific)
    savefig("total_energy_error.pdf")
end
