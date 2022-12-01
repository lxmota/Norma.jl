function loop(params::Dict{Any, Any})
    model = params["model_struct"]
    solver = params["solver_struct"]
    num_steps = params["number of steps"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_diff = final_time - initial_time
    time_step = time_diff / num_steps
    initialize_writing(model)
    for stop ∈ 0 : num_steps
        time = initial_time + stop * time_step
        model.time = time
        apply_bcs(model)
        solve(model, solver)
        write_step(model, stop + 1, time)
        println("Stop ", stop, " of ", num_steps, ", time: ", time)
    end
    finalize_writing(model)
end