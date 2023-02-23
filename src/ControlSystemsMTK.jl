module ControlSystemsMTK
using RobustAndOptimalControl: NamedStateSpace
#=
Ideas: All connections handled by ModelingToolkit.
Names: 
- handled either by named system, or directly in constructor to ODESystem. 
Functions: 
- Give me linear system from [u1, u3] to [qm, a]
If the linearization of a full system produces a named system, one could implement getindex for vectors of names and obtain the desired transfer functions.


Another idea: use modelingtoolkitize/build_function to obtain a function that can be differentiated with ForwardDiff or FD, https://discourse.julialang.org/t/differentialequations-control-systems-and-linearization/31178/6


A third idea: just use named systems with named indexing to obtain any system you want.

=#
using LinearAlgebra
using ModelingToolkit, ControlSystemsBase
using ControlSystemsBase: ssdata, AbstractStateSpace, Continuous, nstates, noutputs, ninputs
# using ControlSystemIdentification
using RobustAndOptimalControl
import ModelingToolkit: ODESystem, FnType, Symbolics
using ModelingToolkit: states, observed, isdifferential
using ModelingToolkit.Symbolics
using ModelingToolkit.Symbolics: jacobian, solve_for
using UnPack
# using Optim, Optim.LineSearches

# using SymbolicControlSystems

export sconnect, feedback, ODESystem, states, observed, named_ss
export batch_ss, trajectory_ss, GainScheduledStateSpace
export build_quadratic_cost_matrix

include("ode_system.jl")
# include("symbolic_optimization.jl")

end
