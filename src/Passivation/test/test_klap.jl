module Test_klap

using Test

using Passivation

using LinearAlgebra, ControlSystemsBase
using FiniteDifferences

@testset "Test_klap.jl" begin
    A = [-1 4; -2 -1];
    B = [1; 2];
    C = [1 0];
    D = 1//8;
    
    Σ = ss(A, B, C, D)

    n, m = Σ.nx, Σ.nu
    M = sqrt(Σ.D + Σ.D')
    P = gram(Σ, :c)

    L = ones(typeof(Σ.A[1,1]), n, m)

    # precompute Schur of A
    F = schur(Σ.A)

    @testset "X" begin
        X = Passivation.X(Σ, L)
        @test norm(Σ.A' * X + X * Σ.A + L * L') < 1e-10
        X = Passivation.X(F, L)
        @test norm(Σ.A' * X + X * Σ.A + L * L') < 1e-10
    end

    @testset "Cr" begin
        Cr = Passivation.Cr(Σ, L, M)
        @test norm(Passivation.X(Σ, L) * Σ.B - Cr' + L*M) < 1e-10
        Cr = Passivation.Cr(Σ, F, L, M)
        @test norm(Passivation.X(F, L) * Σ.B - Cr' + L*M) < 1e-10
    end

    @testset "J" begin
        Σp = Passivation.Σp(Σ, L, M)
        @test norm(Passivation.J(Σ, L, M, P) - norm(Σ - Σp)^2) < 1e-12
    end

    @testset "∇J" begin
        f(L) = Passivation.J(Σ, L, M, P)
        g(L) = Passivation.∇J(Σ, L, M, P)
        @test norm(grad(central_fdm(5,1), f, L)[1] - g(L)) < 1e-8
    end

    @testset "fg!" begin
        f = Inf
        G = zero(L)
        
        f = Passivation.fg!(f, G, L, Σ, P, M)
        @test norm(f - Passivation.J(Σ, L, M, P)) < 1e-10
        @test norm(G - Passivation.∇J(Σ, L, M, P)) < 1e-10

        f = Inf
        G = zero(L)

        f = Passivation.fg!(f, G, L, Σ, F, P, M)
        @test norm(f - Passivation.J(Σ, L, M, P)) < 1e-10
    end

    @testset "klap" begin
        L0 = [-1.5; 1.5;;]
        Σp, _ = klap(Σ, L0)
        @test norm(Σ - Σp)^2 ≈ 2.5

        L0 = [1.5; 1.5;;]
        Σp, _ = klap(Σ, L0)
        @test norm(Σ - Σp)^2 ≈ 0.12751294803962673

        L0 = [-1.5; 1.5;;]
        result = klap(Σ, L0; restart=true)
        @test length(result) == 3
        Σp = result[1]
        @test norm(Σ - Σp)^2 ≈ 0.12751294803962673
    end
end

end
