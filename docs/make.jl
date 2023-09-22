ENV["GKSwstype"] = 322 # workaround for gr segfault on GH actions
# ENV["GKS_WSTYPE"]=100 # try this if above does not work
using Documenter, ControlSystemsMTK, RobustAndOptimalControl, ControlSystemsBase, ModelingToolkit, ModelingToolkitStandardLibrary, MonteCarloMeasurements, SymbolicControlSystems

using Plots
gr()


makedocs(
      sitename = "ControlSystemsMTK Documentation",
      doctest = false,
      modules = [ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, ModelingToolkitStandardLibrary, RobustAndOptimalControl, SymbolicControlSystems],
      remotes = Dict(
            dirname(dirname(pathof(ControlSystemsMTK))) => (Remotes.GitHub("JuliaControl", "ControlSystemsMTK.jl"), "1"),
            dirname(dirname(pathof(ControlSystemsBase))) => (Remotes.GitHub("JuliaControl", "ControlSystems.jl"), "1"),
            dirname(dirname(pathof(ModelingToolkit))) => (Remotes.GitHub("SciML", "ModelingToolkit.jl"), "8"),
            dirname(dirname(pathof(ModelingToolkitStandardLibrary))) => (Remotes.GitHub("SciML", "ModelingToolkitStandardLibrary.jl"), "2"),
            dirname(dirname(pathof(RobustAndOptimalControl))) => (Remotes.GitHub("JuliaControl", "RobustAndOptimalControl.jl"), "0.4"),
            dirname(dirname(pathof(SymbolicControlSystems))) => (Remotes.GitHub("JuliaControl", "SymbolicControlSystems.jl"), "0.1"),
      ),
      pages = [
            "Home" => "index.md",
            "Examples" => [
                  "Batch linearization and gain scheduling" => "batch_linearization.md",
            ],
            "API" => "api.md",
      ],
      format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
      warnonly = [:missing_docs, :docs_block, :cross_references],
)

deploydocs(
      repo = "github.com/JuliaControl/ControlSystemsMTK.jl.git",
)
