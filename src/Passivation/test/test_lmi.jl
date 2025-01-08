module Test_lmi

using Test

using Passivation

using LinearAlgebra, ControlSystemsBase, JuMP

@testset "test_lmi.jl" begin
    A = [-1 4; -2 -1];
    B = [1; 2];
    C = [1 0];
    D = 1//8;
    
    Σ = ss(A, B, C, D)
    
    @testset "lmi" begin
        Σp, model = lmi(Σ)
        X = value.(object_dictionary(model)[:X])
        @test isposdef(X)
        @test isposdef(kypmat(Σp, X) + 1e-6 * I)
        @test norm(Σp - Σ)^2 < 2e-1
    end

    @testset "lmi_tp" begin
        Σp, model = lmi_tp(Σ)
        W = value.(object_dictionary(model)[:W])
        @test isposdef(W + 1e-6 * I)
        @test norm(Σp - Σ)^2 < 2e-1

        Σpref, model = lmi(Σ)
        X = value.(object_dictionary(model)[:X])
        @test norm(Σp.C - Σpref.C) < 1e-6
        @test norm(W - kypmat(Σpref, X)) < 1e-6
    end
    
end

end