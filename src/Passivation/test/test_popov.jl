module Test_popov

using Test

using Passivation

using LinearAlgebra, ControlSystemsBase

@testset "test_popov.jl" begin
    A = [-2 1.; -1. -1.]
    B = [6.; 0]
    C = B'
    D = 1
    Σp = ss(A,B,C,D)

    A = [-1 4; -2 -1];
    B = [1; 2];
    C = [1 0];
    D = 1//8;
    Σnp = ss(A, B, C, D)

    @testset "popov" begin
        @testset "s" begin
            s = 1 + 1im
            @test norm(popov(Σp, s) - (evalfr(Σp, s) + transpose(evalfr(Σp, -s)))) < 1e-8
        end

        @testset "ω" begin
            ω = 10 .^ range(-2, stop=2, length=5)
            @test norm(popov(Σp, ω) - (freqresp(Σp, ω) + permutedims(freqresp(Σp, -ω), [2,1,3]))) < 1e-8
        end            
    end

    @testset "test_ispassive" begin
        for opt ∈ [:lmi, :popov]
            @test ispassive(Σp; opt=opt) == true
            @test ispassive(Σnp; opt=opt) == false
        end
    end
end

end