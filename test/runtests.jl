using ControlSystemsMTK
using Test

@testset "ControlSystemsMTK.jl" begin
    @testset "ODESystem" begin
        @info "Testing ODESystem"        
        include("test_ODESystem.jl")
    end

    @testset "batchlin" begin
        @info "Testing batchlin"
        include("test_batchlin.jl")
    end

    # @testset "symbolic_opt" begin
    #     @info "Testing symbolic_opt"
    #     include("test_symbolic_opt.jl")
    # end

end
