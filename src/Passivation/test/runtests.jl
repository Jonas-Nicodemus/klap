using Test

@testset "Passivation.jl" begin
    include("test_lyapunov.jl")
    include("test_klap.jl")
    include("test_lmi.jl")
    include("test_kyp.jl")
    include("test_popov.jl")
end