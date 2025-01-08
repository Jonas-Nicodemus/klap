"""
    X = kyp(Σ; optimizer=Hypatia.Optimizer, kwargs...)

Returns a solution to the KYP inequality by solving the corresponding linear matrix inequality.
"""
function kyp(Σ::StateSpace; kwargs...)
    model = Model(() -> Hypatia.Optimizer(verbose = false))
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, X[1:Σ.nx, 1:Σ.nx], PSD)
    @constraint(model, kypmat(Σ, X) in PSDCone())
    optimize!(model)

    return value.(X)
end

"""
    Xmin = kypmin(Σ; kwargs...)

Returns the minimal solution to the KYP inequality by solving the Riccati equation for the stabilizing solution.
"""
function kypmin(Σ; kwargs...)
    return kypextremal(Σ, true, kwargs...)
end

"""
    Xmin = kypmax(Σ; kwargs...)

Returns the maximal solution to the KYP inequality by solving the Riccati equation for the anti-stabilizing solution.
"""
function kypmax(Σ; kwargs...)
    return kypextremal(Σ, false, kwargs...)
end

"""
    X = kypmax(Σ; kwargs...)

Returns the maximal or minimal solution to the KYP inequality by solving the Riccati equation for the stabilizing or anti-stabilizing solution.
"""
function kypextremal(Σ, min=true; kwargs...)
    A = Array(Σ.A)
    B = Array(Σ.B)
    R = Array(-Σ.D - Σ.D')
    Q = zero(Σ.A)
    S = Array(-Σ.C')

    X, _, _ = arec(A, B, R, Q, S, as=!min)
    
    return X
end

"""
    W = kypare(Σ, X)
    
Returns the residual of the KYP ARE for some solution candidate `X`.
"""
function kypare(Σ::StateSpace, X)
   return Σ.A' * X + X * Σ.A + (Σ.C' - X * Σ.B) * inv(Σ.D + Σ.D') * (Σ.C - Σ.B' * X)
end

"""
    W = kypmat(Σ, X) 

Returns the KYP matrix of the system `Σ` for the given matrix `X`.
"""
function kypmat(Σ::StateSpace, X)
    return [-Σ.A'*X-X*Σ.A  Σ.C'-X*Σ.B; Σ.C-Σ.B'*X Σ.D + Σ.D']
end