module Passivation

using LinearAlgebra, ControlSystemsBase, MatrixEquations, Plots
using Optim, JuMP, Hypatia

export passivate, klap, klap_inital_guess, lmi, lmi_tp
export ispassive, popov, popovplot, sampopov
export kyp, kypmin, kypmax, kypare, kypmat

include("lyapunov.jl")
include("klap.jl")
include("lmi.jl")
include("kyp.jl")
include("popov.jl")

"""
    passivate(Σ::StateSpace, method=:klap, args...; kwargs...) -> Σp, res

Passivate the system `Σ` using the method `method`. The available methods are:

- `:klap`: KLAP optimization
- `:lmi`: LMI optimization [GS21](@cite)
- `:lmi_tp`: LMI optimization with trace parametrization [Dum02, CPS04](@cite)

The remaining arguments `args` and keyword arguments `kwargs` are passed to the corresponding passivation function.
"""
function passivate(Σ::StateSpace, method=:klap, args...; kwargs...)
    if method == :klap
        return klap(Σ, args...; kwargs...)
    elseif method == :lmi
        return lmi(Σ, args...; kwargs...)
    elseif method == :lmi_tp
        return lmi_tp(Σ, args...; kwargs...)
    else
        error("Unknown passivation method $method")
    end
end

end