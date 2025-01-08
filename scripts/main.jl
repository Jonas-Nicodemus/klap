using DrWatson
@quickactivate "KLAP"

using LinearAlgebra, ControlSystemsBase, RobustAndOptimalControl
using DataFrames, Plots, LaTeXStrings, JuMP, Optim
using Passivation

import Passivation: ispassive

# ACC or CDplayer
model_name = "ACC"; # "ACC", "CDplayer"
model_path = datadir(model_name * ".jld2");

Σ = load(model_path)["Σ"];
ispassive(Σ, opt=:popov)
popovplot(Σ)

methods = [:lmi, :lmi_tp, :klap];
results = DataFrame(method=Symbol[], solve_time=Float64[], iterations=Int[], h2_error=Float64[])

for method in methods
    Σp, res = passivate(Σ, method; verbose=false)
    if isa(res, JuMP.Model)
        res = solution_summary(res);
        solve_time = res.solve_time
        iterations = res.barrier_iterations
        # h2_error = res.objective_value
        
    elseif isa(res, Optim.OptimizationResults)
        solve_time = res.time_run
        iterations = res.iterations
        # h2_error = sqrt(res.minimum)
    else
        error("Unknown result type")
    end

    h2_error = h2norm(Σ - Σp)
    
    @info "Method: $method,\t solve time: $solve_time,\t iterations: $iterations,\t H2-error: $h2_error"
    push!(results, (method, solve_time, iterations, h2_error))
end

# print results (Note: Due to Julia’s nature, rerun from line 19 to obtain results excluding compilation time.)
@info results
