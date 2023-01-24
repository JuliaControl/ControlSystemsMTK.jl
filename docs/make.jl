ENV["GKSwstype"] = 322 # workaround for gr segfault on GH actions
# ENV["GKS_WSTYPE"]=100 # try this if above does not work
using Documenter, ControlSystemsMTK, RobustAndOptimalControl, ControlSystemsBase, ModelingToolkit, ModelingToolkitStandardLibrary, MonteCarloMeasurements, SymbolicControlSystems

using Plots
gr()


makedocs(
      sitename = "ControlSystemsMTK Documentation",
      doctest = false,
      modules = [ControlSystemsMTK, ControlSystemsBase, ModelingToolkit, ModelingToolkitStandardLibrary, RobustAndOptimalControl, SymbolicControlSystems],
      strict=[
        :doctest, 
        :linkcheck, 
        :parse_error,
        :example_block,
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
      ],
      pages = [
            "Home" => "index.md",
            "Examples" => [
                  "Batch linearization and gain scheduling" => "batch_linearization.md",
            ],
            "API" => "api.md",
      ],
      format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
)

deploydocs(
      repo = "github.com/JuliaControl/ControlSystemsMTK.jl.git",
)
