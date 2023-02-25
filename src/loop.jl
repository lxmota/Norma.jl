function loop(params::Dict{Any,Any})
    sim_type = params["type"]
    if sim_type == "single"
        return loop_single(params)
    elseif sim_type == "multi"
        return loop_multi(params)
    else
        error("Unknown type of simulation: ", sim_type)
    end
end

function loop_single(params::Dict{Any,Any})
    model = create_model(params)
    solver = create_solver(params)
    integrator = create_time_integrator(params)
    apply_ics(model)
    initialize_writing(integrator, model)
    while integrator.time <= integrator.final_time
        model.time = integrator.time
        println("Stop ", integrator.stop, ", time: ", integrator.time)
        apply_bcs(model)
        update_dofs(model, solver)
        solve(integrator, model, solver)
        write_step(integrator, model)
        integrator.time = round(integrator.time + integrator.time_step; digits=10)
        integrator.stop += 1
    end
    finalize_writing(model)
    return integrator, solver, model
end

function loop_multi(params::Dict{Any,Any})
    models = Vector{Model}()
    integrators = Vector{TimeIntegrator}()
    solvers = Vector{Solver}()
    domains = params["domains"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    for domain ∈ domains
        domain_params = params[domain]
        integrator_params = domain_params["time integrator"]
        integrator_params["initial time"] = initial_time
        integrator_params["final time"] = final_time
        model = create_model(domain_params)
        solver = create_solver(domain_params)
        integrator = create_time_integrator(domain_params)
        push!(models, model)
        push!(solvers, solver)
        push!(integrators, integrator)
    end
    return integrators, solvers, models
end