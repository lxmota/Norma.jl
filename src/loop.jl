function update_dofs(model::Any, solver::Any)
    solver.free_dofs = model.nodal_dofs .== free::DOF
end

function loop(params::Dict{Any, Any})
    integrator = params["time_integrator_struct"]
    model = params["model_struct"]
    solver = params["solver_struct"]
    apply_ics(model)
    initialize_writing(model)
    while integrator.time <= integrator.final_time
        model.time = integrator.time
        println("Stop ", integrator.stop, ", time: ", integrator.time)
        apply_bcs(model)
        update_dofs(model, solver)
        solve(integrator, model, solver)
        write_step(model, integrator.stop + 1, integrator.time)
        integrator.time += integrator.time_step
        integrator.stop += 1
    end
    finalize_writing(model)
    return model, solver
end