function StaticSchwarzController(params::Dict{Any,Any})
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    time = initial_time
    stop =0
    previous_solutions = Vector{Vector{Float64}}()
    StaticSchwarzController(initial_time, final_time, time_step, time, stop, previous_solutions)
end

function DynamicSchwarzController(params::Dict{Any,Any})
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    time = initial_time
    stop =0
    previous_solutions = Vector{Vector{Float64}}()
    previous_sol_dots = Vector{Vector{Float64}}()
    previous_sol_dotdots = Vector{Vector{Float64}}()
    DynamicSchwarzController(initial_time, final_time, time_step, time, stop, previous_solutions, previous_sol_dots, previous_sol_dotdots)
end

function create_schwarz_controller(params::Dict{Any,Any})
    type = params["subdomains type"]
    if type == "static"
        return StaticSchwarzController(params)
    elseif type == "dynamic"
        return DynamicSchwarzController(params)
    else
        error("Unknown type of Schwarz controller : ", type)
    end
end

function schwarz(simulation::MultiDomainSimulation, initial_time::Float64, final_time::Float64, stop::Int64)
    apply_bcs(simulation)
    if stop == 0
        solve(simulation)
        write_step(simulation)
    end
    write_step(simulation)
end