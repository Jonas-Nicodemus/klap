"""
    lmi_tp(Σ::StateSpace; kwargs...) -> Σp, model

Solves the passivation problem for a given state-space system `Σ` using the positive real LMI constraints with trace parametrization [Dum02, CPS04](@cite).
"""
function lmi_tp(Σ::StateSpace; kwargs...)
    # Dimensions
    n = Σ.nx
    m = Σ.nu

    Qc = cholesky(gram(Σ, :c)).U

    ζ = zeros(n, n, n, m)
    for i in 1:n
        for j in 1:m
            ζ[:,:,i,j] = lyapc(Σ.A, hermitianpart(Σ.B[:,j] * I[1:n, i]'))
        end
    end

    # Passivity Enforcement
    model = Model(Hypatia.Optimizer)
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, W[1:n+m, 1:n+m], PSD)
    fix.(W[n+1:end, n+1:end], Σ.D + Σ.D')

    Ψ1 = W[1:n, 1:n]
    Ψ2 = W[1:n, n+1:end]

    @variable(model, t)
    @variable(model, Δ[1:m, 1:n])
    @constraint(model, c[i=1:n, j=1:m], Σ.C[j,i] + Δ[j,i] == tr(ζ[:,:,i,j] * Ψ1) + Ψ2[i, j])
    @constraint(model, [t; vec(Δ * Array(Qc'))] in SecondOrderCone())
    @objective(model, Min, t)

    JuMP.optimize!(model)

    return ss(Σ.A, Σ.B, Σ.C + value.(Δ), Σ.D), model;
end

"""
    lmi(Σ::StateSpace; kwargs...) -> Σp, model

Solves the passivation problem for a given state-space system `Σ` using the positive real LMI constraints [GS21](@cite).
"""
function lmi(Σ::StateSpace; kwargs...)
    # Dimensions
    n = Σ.nx
    m = Σ.nu
    
    # Passivity Enforcement
    Wc = gram(Σ, :c)
    Qc = cholesky(Wc).U
    
    model = Model(Hypatia.Optimizer)
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, X[1:n, 1:n], PSD)
    @variable(model, Δ[1:m, 1:n])

    C = Σ.C + Δ
    W = -[Σ.A'*X+X'*Σ.A X*Σ.B-C'; Σ.B'*X-C -Σ.D-Σ.D']

    @variable(model, t)
    @objective(model, Min, t)
    @constraint(model, [t; vec(Δ * Array(Qc'))] in SecondOrderCone())
  
    @constraint(model, W in PSDCone())
    
    JuMP.optimize!(model)
    
    return ss(Σ.A, Σ.B, value.(C), Σ.D), model
end