function loop(params::Dict{Any, Any})
    num_steps = params["number of steps"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_diff = final_time - initial_time
    for stop ∈ 0 : num_steps
        time = initial_time + stop * time_diff 
    end
end