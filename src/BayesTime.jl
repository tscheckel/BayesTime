module BayesTime
    # Import or export symbols
    export forcePD
    export get_companion
    export get_hs
    export get_IRF
    export mlag
    export offsetting
    export CCCM
    export SAVS
    export trace_plot
    export lower_tri
    export VAR_OLS
    export diff_lag
    export data_trans
    
    # globally required packages
    using Distributions
    using LinearAlgebra
    using MCMCDiagnosticTools
    using Plots
    using DataFrames

    # Include submodules or helper files
    include("get_hs.jl")
    include("mlag.jl")
    include("offsetting.jl")
    include("get_companion.jl")
    include("get_IRF.jl") # requires get_companion.jl
    include("forcePD.jl")
    include("CCCM.jl")
    include("SAVS.jl")
    include("trace_plot.jl")
    include("lower_tri.jl")
    include("VAR_OLS.jl")
    include("diff_lag.jl")
    include("data_trans.jl") # requires diff_lag.jl
end