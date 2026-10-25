"""
    get_VAR_design(; data, Exo = nothing, P = 1, cons = true, trend = false)

Create the dependent-variable and design matrices of a VAR(P).

The design matrix is X = [Y_{t-1} ... Y_{t-P} Exo trend const], where exogenous variables,
trend and constant are appended only if requested. The first P observations are lost to the lags.

# Arguments
- `data`: TxM data matrix or `DataFrame` (rows are observations)
- `Exo`: optional TxN matrix or `DataFrame` of exogenous variables (default: nothing)
- `P::Int`: number of lags (default: 1)
- `cons::Bool`: include a constant (default: true)
- `trend::Bool`: include a linear time trend (default: false)

# Returns
- `Dict` containing:
    - `:Y`: (T-P)xM matrix of dependent variables
    - `:X`: (T-P)xK design matrix
    - `:reg_names`: names of the K regressors
    - `:var_names`: names of the M variables

# Author
Tobias Scheckel
"""
function get_VAR_design(;
    data::Union{AbstractMatrix{<:Real}, DataFrame}, Exo::Union{AbstractMatrix{<:Real}, DataFrame, Nothing} = nothing,
    P::Int = 1, cons::Bool = true, trend::Bool = false
)
    # CHECK INPUT DATATYPES
    if isa(data, DataFrame)
        data = NamedArray(
            Matrix{Float64}(data),
            (collect(axes(data, 1)), collect(names(data))),
            ("obs", "vars")
        )
    end
    if isa(Exo, DataFrame)
        Exo = NamedArray(
            Matrix{Float64}(Exo),
            (collect(axes(Exo, 1)), collect(names(Exo))),
            ("obs", "exo")
        )
    end

    var_names = collect(1:size(data, 2))
    if data isa NamedArray
        var_names = collect(names(data, 2))
    end

    # create matrix of dependent variables and regressors
    Y = Matrix{Float64}(data[P+1:end,:])
    T = size(Y, 1)
    # create design matrix
    X_lag = mlag(data, P)
    X = Matrix{Float64}(X_lag)
    reg_names = collect(1:size(X, 2))
    if X_lag isa NamedArray
        reg_names = collect(names(X_lag, 2))
    end
    # add exogenous regressors if provided
    if !isnothing(Exo)
        Exo_slice = Exo[P+1:end,:]
        X = hcat(X, Matrix{Float64}(Exo_slice))
        exo_names = ["exo_$i" for i in 1:size(Exo_slice, 2)]
        if Exo_slice isa NamedArray
            exo_names = collect(names(Exo_slice, 2))
        end
        reg_names = vcat(reg_names, exo_names)
    end
    if trend
        X = hcat(X, collect(1.0:T))
        reg_names = vcat(reg_names, ["trend"])
    end
    if cons
        X = hcat(X, fill(1.0, T))
        reg_names = vcat(reg_names, ["const"])
    end

    retDict = Dict(
        :Y => Y,
        :X => X,
        :reg_names => reg_names,
        :var_names => var_names
    )

    return retDict
end
