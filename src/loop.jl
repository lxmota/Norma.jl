function update_dofs(model::Any, solver::Any)
    solver.free_dofs = model.nodal_dofs .== free::DOF
end

function loop(params::Dict{Any, Any})
    time_integrator = params["time_integrator_struct"]
    model = params["model_struct"]
    solver = params["solver_struct"]
    time = time_integrator.initial_time
    stop = 0
    initialize_writing(model)
    while time <= time_integrator.final_time
        model.time = time
        println("Stop ", stop, ", time: ", time)
        apply_bcs(model)
        update_dofs(model, solver)
        solve(time_integrator, model, solver)
        write_step(model, stop + 1, time)
        time += time_integrator.time_step
        stop += 1
    end
    finalize_writing(model)
    return model, solver
end