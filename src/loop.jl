include("apply_bcs.jl")

function loop(params::Dict{Any, Any})
    model = params["model_struct"]
    mesh_struct = params["mesh_struct"]
    num_steps = params["number of steps"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_diff = final_time - initial_time
    time_step = time_diff / num_steps
    for stop ∈ 0 : num_steps
        time = initial_time + stop * time_step
        model.time = time
        apply_bcs(model)
        energy = potential_energy(model)
        println("Time: ", time, ", Energy: ", energy)
    end
end