module BayesTime
    # Import or export symbols
    export forcePD
    export get_companion
    export get_hs
    export get_IRF
    export mlag
    export offsetting
    export CCCM
    
    # globally required packages
    using Distributions
    using LinearAlgebra

    # Include submodules or helper files
    include("get_hs.jl")
    include("mlag.jl")
    include("offsetting.jl")
    include("get_companion.jl")
    include("get_IRF.jl") # requires get_companion.jl
    include("forcePD.jl")
    include("CCCM.jl")
end # module