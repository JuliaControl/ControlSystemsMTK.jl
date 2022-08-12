using ControlSystemsMTK
using Test

@testset "ControlSystemsMTK.jl" begin
    @testset "ODESystem" begin
        @info "Testing ODESystem"        
        include("test_ODESystem.jl")
    end

    # @testset "symbolic_opt" begin
    #     @info "Testing symbolic_opt"
    #     include("test_symbolic_opt.jl")
    # end

end
