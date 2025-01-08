"""
    ispassive(Σ; opt=:lmi, kwargs...)

Checks if the system `Σ` is passive. The function returns `true` if the system is passive and `false` otherwise.
If `opt` is set to `:lmi`, the KYP lmi is tried to solve. If `opt` is set to `:popov`, the Popov function is sampled to check passivity.
"""
function ispassive(Σ::StateSpace; opt=:lmi, kwargs...)
    if opt==:lmi
        return ispassive_lmi(Σ; kwargs...)
    elseif opt==:popov
        return ispassvie_sampopov(Σ; kwargs...)
    else 
        error("opt must be either :lmi or :popov")
    end
end

"""
    popov(Σ, s)
    popov(Σ, ω)

Evaluates the popov function of the system `Σ` at the complex variable `s` or real frequencies `ω`.

    Φ(s) = Σ(s) + Σ(-s)'

where `Σ(s)` is the frequency response (transfer function) of `Σ` at the complex variable `s`.
"""
function popov(Σ::StateSpace, s::Complex)
    return evalfr(Σ, s) + transpose(evalfr(Σ, -s))
end

function popov(Σ::StateSpace, ω::AbstractVector{W}) where W <: Real
    return popov(freqresp(Σ, ω))
end

function popov(H)
    return real(H + permutedims(conj(H), [2,1,3]))
end

"""
    sampopov(Σ; ω=10 .^ range(-15, stop=5, length=5000)) -> Λmin, ω, Ωp, Ωnp, F, c

Samples the Popov function of the system `Σ` at the frequencies `ω`.
Returns the minimum eigenvalues of the sampled Popov function, the frequencies `ω`, the passive and non-passive intervals `Ωp` and `Ωnp`, the eigendecompositions `F` and the change points `c`.
"""
function sampopov(Σ; ω=10 .^ range(-15, stop=5, length=5000))
    Φ = popov(Σ, ω)
    F = eigen.([Φ[:, :, i] for i ∈ axes(Φ, 3)])
    Λmin = [minimum(real.(F[i].values)) for i ∈ eachindex(F)]
    
    c = findall(i -> (sign(Λmin[i]) != sign(Λmin[i+1])), range(1, length(Λmin)-1))
    c = [0; c...; length(Λmin)]
    Ωi = [c[i]+1:c[i+1] for i ∈ range(1, length(c)-1)]
    
    if Λmin[1] > 0
        Ωp = Ωi[1:2:end]
        Ωnp = Ωi[2:2:end] 
    else
        Ωp = Ωi[2:2:end]
        Ωnp = Ωi[1:2:end]
    end

    return Λmin, ω, Ωp, Ωnp, F, c
end

function ispassvie_sampopov(Σ; kwargs...)
    _, _, _, Ωnp, _, _  = sampopov(Σ; kwargs...)
    return isempty(Ωnp)
end

function ispassive_lmi(Σ::StateSpace; ε=1e-8, kwargs...)
    X = kyp(Σ; kwargs...)
    return isposdef(X) && isposdef(kypmat(Σ, X) + ε * I)
end

"""
    popovplot(Σ; args...; kwargs...)
    popovplot!(p, Σ; args...; kwargs...)

Plots the minimal eigenvalue of the Popov function of the system `Σ` at the frequencies `ω`.
"""
function popovplot(Σ, args...; kwargs...)
    p = plot()
    
    popovplot!(p, Σ, args...; kwargs...)

    return p
end

function popovplot!(Σ, args...; kwargs...)
    popovplot!(current(), Σ, args...; kwargs...)
end

function popovplot!(p, Σ, args...; kwargs...)
    Λmin, ω, Ωp, Ωnp, _, _ = sampopov(Σ; kwargs...)

    plot!(p, label="min(λ)", xaxis=:log, legend=nothing,
        xlabel="Frequency ω", ylabel="Φ(ω)", title="Popov plot")

    hline!(p, [0], color=:black, linestyle=:dash, label=nothing)

    for Ωi ∈ Ωp
        plot!(p, ω[Ωi], Λmin[Ωi]; color=:green, label=nothing)
    end

    for Ωi ∈ Ωnp
        plot!(p, ω[Ωi], Λmin[Ωi]; color=:red, label=nothing)
    end
end