module BayesTime
    # Import or export symbols
    export get_hs

    # Include submodules or helper files
    include("get_hs.jl")
    include("mlag.jl")
    include("offsetting.jl")
    include("get_companion.jl")
    include("get_IRF.jl") # requires get_companion.jl

end # module