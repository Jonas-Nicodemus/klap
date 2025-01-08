using DrWatson, Test
@quickactivate "klap"

# Run test suite
println("Starting tests")
ti = time()

@testset "KLAP" begin
    include(srcdir("Passivation", "test", "runtests.jl"))
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")