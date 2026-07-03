module ControlSystemsMTK
using RobustAndOptimalControl: NamedStateSpace
#=
Ideas: All connections handled by ModelingToolkit.
Names: 
- handled either by named system, or directly in constructor to System. 
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
using RobustAndOptimalControl, MonteCarloMeasurements
import ModelingToolkit: System, FnType, Symbolics
using ModelingToolkit: unknowns, observed, isdifferential
using Symbolics
using Symbolics: jacobian, solve_for
using UnPack
# using Optim, Optim.LineSearches

# using SymbolicControlSystems

export feedback, System, unknowns, observed, named_ss
export batch_ss, trajectory_ss, GainScheduledStateSpace
export build_quadratic_cost_matrix

export get_named_sensitivity, get_named_comp_sensitivity, get_named_looptransfer

"""
Default keyword arguments forwarded to ModelingToolkit's linearization (and through it to
`mtkcompile`) by every linearization-based function in this package (`named_ss`,
`get_named_sensitivity` and friends, `batch_ss`, `trajectory_ss`,
`build_quadratic_cost_matrix`). `inline_linear_sccs` makes tearing solve linear
strongly-connected components inline (symbolically for size ≤ `analytical_linear_scc_limit`)
rather than through numerical linear-solve calls embedded in the generated code. This is
required for correct linearization of models with large linear SCCs (e.g. multibody
mechanisms). Keyword arguments passed by the caller take precedence, so pass
`reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm()` to
restore the ModelingToolkit default.
"""
const DEFAULT_LINEARIZE_KWARGS = (;
    reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(;
        inline_linear_sccs = true, analytical_linear_scc_limit = 1),
)

include("ode_system.jl")
# include("symbolic_optimization.jl")

end
