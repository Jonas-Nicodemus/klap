"""
    klap(Σ::StateSpace; L0=L0(Σ), M=M(Σ), P=gram(Σ, :c); recycl=:schur, restart=false, α=1e-8, ε=1e-4, verbose=true, kwargs...) -> Σp, res

Passivates a system `Σ` using the KLAP method. The optimization problem is solved using LBFGS.
"""
function klap(Σ::StateSpace, L0=L0(Σ), M=M(Σ), P=gram(Σ, :c); recycl=:schur, restart=false, α=1e-8, ε=1e-4, verbose=true, kwargs...)
    if isdiag(Σ.A)
        @info "Diagonal Σ.A detected, using specialized fg! implementation"
        Ã = diag(Σ.A) .+ diag(Σ.A)'
        d = Optim.only_fg!((f,G,L) -> fg!(f, G, L, Σ, Diagonal(Σ.A), Ã, P, M))
    elseif recycl == :schur
        S = schur(Σ.A)
        d = Optim.only_fg!((f,G,L) -> fg!(f, G, L, Σ, S, P, M))
    elseif recycl === nothing
        @info "No recycling method specified, using generic fg! implementation"
        d = Optim.only_fg!((f,G,L) -> fg!(f, G, L, Σ, P, M))
    else 
        @error "Recycling method $recycl not recognized"
    end
    
    res = Optim.optimize(d, L0, LBFGS(),  
        Optim.Options(;
            g_abstol = 1e-8,
            g_reltol = 1e-8,
            f_tol = 1e-6,
            x_tol = 0.0,
            show_trace = verbose,
            show_every = 1,
            extended_trace = false,
            store_trace = false,
            iterations = 10_000,
            kwargs...
        ))

    if verbose
        @info "Optimization result: $(res)"
    end

    L = Optim.minimizer(res)
    Σp = Passivation.Σp(Σ, L, M)

    if restart 
        Y = Passivation.Y(Σp, L, M)
        @info "Checking eigvals of Y ..."
       
        if norm(real(eigvals(Y))) / norm(Σ.A) > ε
            @info "Possible stuck in local minimum, performing unconstrained gradient step ..."
           
            Gc = 2*(Σp.C-Σ.C)*P
            Σtmp = ss(Σ.A, Σ.B, Σp.C - α*Gc, Σ.D)
            L, ΔD = Passivation.klap_inital_guess(Σtmp)

            if ΔD > 0
                @warn "The system was not passive after the unconstrained gradient step. Perturbation ΔD = $ΔD was used.
                Consider decreasing α or increasing ε."
            end

            return klap(Σ, L, M, P; recycl=recycl, restart=restart, α=α, verbose=verbose, kwargs...)..., res
        else
            @info "Global minimum detected, no restart needed."
        end
    end

    return Σp, res
end

"""
    klap_inital_guess(Σ, ΔD=0.0; ε=0) -> L0, ΔD

Computes an initial guess for the KLAP optimization problem. The initial guess is computed by perturbing the feedthrough matrix to achieve a passive realization.
Then the perturbed system is used to compute the initial guess. The perturbation `ΔD` can be specified, otherwise it is computed using `ΔD(Σ)`.
"""
function klap_inital_guess(Σ, ΔD=0.0; ε=0)
    ΔD = ΔD > 0 ? ΔD : Passivation.ΔD(Σ, ε=ε)

    Σpert = ss(Σ.A, Σ.B, Σ.C, Σ.D + ΔD * I)
    return Passivation.L(Σpert, Passivation.M(Σpert)), ΔD
end

function ΔD(Σ; ε=0)
    Λmin, _, _, _, _, _ = sampopov(Σ)
    ΔD = maximum([0, -minimum(Λmin)/2]) + ε    
    return ΔD
end

function L0(Σ; ε=0)
    return klap_inital_guess(Σ, ε=ε)[1]
end

function L(Σ, M=M(Σ))
    Xmin = kypmin(Σ)
    return (Σ.C' - Xmin * Σ.B)*inv(M)
end

function M(Σ)
    return sqrt(Σ.D + Σ.D')
end

function J(Σ::StateSpace, L, M=sqrt(Σ.D + Σ.D'), P=gram(Σ, :c))
    Ce = Σ.C - Cr(Σ, L, M)
    return J(Ce, P)
end

function J(Ce, P)
    return real(tr(Ce * (P * Ce')))
end

function X(F::LinearAlgebra.Factorization, L)
    return lyapc(F, L*L'; adj=true)
end

function X(Ã, U, U⁻¹, L)
    return lyapc(Ã, U, U⁻¹, L*L'; adj=true)
end

function X(A::LinearAlgebra.Diagonal, Ã, L)
    return lyapc(A, L*L', Ã; adj=true)
end

function X(Σ, L)
    U = plyapc(Array(Σ.A'), L)
    return U * U'
end

function Y(Σ, X)
    return Σ.A - Σ.B * inv(Σ.D + Σ.D') * (Σ.C - Σ.B' * X)
end

function Y(Σ, L, M)
    return Σ.A - Σ.B * inv(Σ.D + Σ.D') * M * L'
end

function Cr(Σ, L, M=sqrt(Σ.D + Σ.D'))
    return Σ.B' * X(Σ, L) + M * L'
end

function Cr(Σ, F::LinearAlgebra.Factorization, L, M=sqrt(Σ.D + Σ.D'))
    return Σ.B' * X(F, L) + M * L'
end

function Cr(Σ, Ã, U, U⁻¹, L, M=sqrt(Σ.D + Σ.D'))
    return Σ.B' * X(Ã, U, U⁻¹, L) + M * L'
end

function Cr(Σ, A::LinearAlgebra.Diagonal, F, L, M=sqrt(Σ.D + Σ.D'))
    return Σ.B' * X(A, F, L) + M * L'
end

function Σp(Σ, L, M=sqrt(Σ.D + Σ.D'))
    return ss(Σ.A, Σ.B, Cr(Σ, L, M), Σ.D)
end

function ∇J(Σ, L, M=sqrt(Σ.D + Σ.D'), P=gram(Σ, :c))
    Cr = Passivation.Cr(Σ, L, M) 
    Ce = Cr - Σ.C
    XX = lyapc(Σ.A, 2*hermitianpart(Σ.B*Ce*P)) 
    return 2XX*L + 2*P*Ce'*M
end

function fg!(f, G, L, Σ, P, M)
    Ce = Cr(Σ, L, M) - Σ.C

    if G !== nothing
        XX = lyapc(Σ.A, 2*hermitianpart(Σ.B*Ce*P))
        G .= 2XX*L + 2*P*Ce'*M
    end
    if f !== nothing
      return J(Ce, P)
    end
end

function fg!(f, G, L, Σ, A::LinearAlgebra.Diagonal, F, P, M)
    Ce = Cr(Σ, A, F, L, M) - Σ.C

    if G !== nothing
        XX = lyapc(A, 2*hermitianpart(Σ.B*Ce*P), F)
        G .= 2XX*L + 2*P*Ce'*M
    end
    if f !== nothing
      return J(Ce, P)
    end
end

function fg!(f, G, L, Σ, F::LinearAlgebra.Factorization, P, M)
    Ce = Cr(Σ, F, L, M) - Σ.C

    if G !== nothing
        XX = lyapc(F, 2*hermitianpart(Σ.B*Ce*P))
        G .= 2XX*L + 2*P*Ce'*M
    end
    if f !== nothing
      return J(Ce, P)
    end
end