module Test_kyp

using Test

using Passivation

using LinearAlgebra, ControlSystemsBase

@testset "test_kyp.jl" begin
    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    A = (J-R)*Q
    B = [6.; 0.;;]
    C = B'
    D = [1.;;]
    
    Σ = ss(A, B, C, D)
        
    @testset "KYP are" begin
        Xmin = kypmin(Σ)
        @test norm(kypare(Σ, Xmin)) < 1e-8

        Xmax = kypmax(Σ)
        @test norm(kypare(Σ, Xmax)) < 1e-8
    end

    @testset "KYP lmi" begin
        X = kyp(Σ)

        @test isposdef(X) 
        @test isposdef(kypmat(Σ, X) + 1e-8 * I)
    end
end

end