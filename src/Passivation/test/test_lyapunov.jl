module Test_lyapunov

using Test

using Passivation

using LinearAlgebra, MatrixEquations, ControlSystemsBase

@testset "Test_lyapunov.jl" begin
    A = [-1 4; -2 -1];
    B = [1; 2];
    C = [1 0];
    D = 1//8;
    Σ = ss(A, B, C, D)
   
    @testset "Schur" begin
        Xref = lyapc(Σ.A, Σ.B * Σ.B')
        Yref = lyapc(Σ.A', Σ.C' * Σ.C)

        S = schur(Σ.A)
        X = lyapc(S, Σ.B * Σ.B')
        @test norm(X - Xref) < 1e-15
        
        Y = lyapc(S, Σ.C' * Σ.C; adj=true)
        @test norm(Y - Yref) < 1e-15
    end

    @testset "Diagonal" begin
        F = eigen(Σ.A)
        # similariry_transform ?
        Σdiag = ss(diagm(F.values), F.vectors \ Σ.B, Σ.C * F.vectors, Σ.D)
        Xref = lyapc(Σdiag.A, Σdiag.B * Σdiag.B')
        Ã = diag(Σdiag.A) * ones(1, size(Σdiag.A,1)) + ones(size(Σdiag.A,1),1) * diag(Σdiag.A)'
        X = lyapc(Diagonal(Σdiag.A), Σdiag.B * Σdiag.B', Ã)
        @test norm(X - Xref) < 1e-15

        Yref = lyapc(Σdiag.A', Σdiag.C' * Σdiag.C)
        Y = lyapc(Diagonal(Σdiag.A), Σdiag.C' * Σdiag.C, Ã; adj=true)
        @test norm(Y - Yref) < 1e-15
    end
    
    end
end