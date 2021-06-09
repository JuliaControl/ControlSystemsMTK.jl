using ControlSystemsMTK
using Test

@testset "ControlSystemsMTK.jl" begin
    @testset "ODESystem" begin
        @info "Testing ODESystem"
        
        include("test_ODESystem.jl")
    end
end
