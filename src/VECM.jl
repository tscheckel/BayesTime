"""
    VECM(; data, Exo=nothing, P::Int=1, r::Int, cons::Bool=true, trend::Bool=false,
          estim::String="ML", map_to_var::Bool=false)

Wrapper around `tsDyn::VECM` via `RCall`.

# Inputs
- `data`: T x M matrix-like object (`AbstractMatrix` or `DataFrame`)
- `Exo`: optional T x K exogenous matrix-like object (`nothing` by default)
- `P`: lag order argument passed to `tsDyn::VECM`
- `r`: cointegration rank
- `cons`: include intercept
- `trend`: include deterministic trend
- `estim`: estimation method in `tsDyn` (default: `"ML"`)
- `map_to_var`: if `true`, also map VECM coefficients into VAR coefficients

# Returns
`Dict` with keys:
- `:B` (`NamedArray`): coefficient matrix from `coef(vecm)`
- `:Pi` (`NamedArray`): long-run impact matrix from `coefPI(vecm)`
- `:cn`, `:rn`: coefficient/equation names from R
- `:A_VAR` (`NamedArray`, only if `map_to_var=true`): mapped VAR coefficients
"""
function VECM(; 
    data::Union{AbstractMatrix{<:Real}, DataFrame},
    Exo::Union{AbstractMatrix{<:Real}, DataFrame, Nothing} = nothing,
    P::Int = 1,
    r::Int,
    cons::Bool = true,
    trend::Bool = false,
    estim::String = "ML",
    map_to_var::Bool = false
)
    P < 1 && throw(ArgumentError("P must be >= 1."))
    r < 0 && throw(ArgumentError("r must be >= 0."))

    data_mat = Matrix{Float64}(data)
    T, M = size(data_mat)

    if !isnothing(Exo)
        exo_mat = Matrix{Float64}(Exo)
        size(exo_mat, 1) == T || throw(ArgumentError("Exo must have the same number of rows as data."))
    else
        exo_mat = nothing
    end

    include = if trend && cons
        "both"
    elseif trend
        "trend"
    elseif cons
        "const"
    else
        "none"
    end

    # Push data into R
    @rput data_mat
    @rput exo_mat
    @rput P
    @rput r
    @rput include
    @rput estim

    R"""
    library(tsDyn)

    if (is.null(exo_mat)) {
      vecm <- VECM(data_mat, lag = P, r = r, include = include, estim = estim)
    } else {
      vecm <- VECM(data_mat, lag = P, r = r, include = include, estim = estim, exogen = exo_mat)
    }

    B <- coef(vecm)
    Pi <- coefPI(vecm)
    cn <- colnames(B)
    rn <- rownames(B)
    """

    # Pull results from R
    @rget B
    @rget Pi
    @rget cn
    @rget rn

    cn = String.(cn)
    rn = String.(rn)

    B_mat = Matrix{Float64}(B)
    Pi_mat = Matrix{Float64}(Pi)

    B_named = NamedArray(B_mat, (rn, cn), ("Eq", "Coeff"))

    # Prefer names from Pi, fall back to data names.
    var_names = if size(Pi_mat, 2) == M
        if Pi isa NamedArray
            String.(names(Pi, 2))
        else
            ["Var" * string(i) for i in 1:M]
        end
    elseif data isa DataFrame
        String.(names(data))
    elseif data isa NamedArray
        String.(names(data, 2))
    else
        ["Var" * string(i) for i in 1:M]
    end

    pi_row_names = rn
    pi_col_names = var_names
    Pi_named = NamedArray(Pi_mat, (pi_row_names, pi_col_names), ("Eq", "LevelVar"))

    ret = Dict{Symbol, Any}(
        :B => B_named,
        :Pi => Pi_named,
        :cn => cn,
        :rn => rn
    )

    if map_to_var
        function gamma_cols_for_lag(lag_idx::Int)
            return [string(v, " -", lag_idx) for v in var_names]
        end

        # Deterministic/exogenous terms that are not ECT and not lagged endogenous terms.
        lag_pattern = Regex(" -[0-9]+$")
        det_exo_cols = [c for c in cn if c != "ECT" && !occursin(lag_pattern, c)]

        A_col_lag = [string(v, "_L", lag_idx) for lag_idx in 1:P for v in var_names]
        A_colnames = vcat(A_col_lag, det_exo_cols)
        A = NamedArray(
            fill(NaN, M, length(A_colnames)),
            (rn, A_colnames),
            ("Eq", "Coeff")
        )

        I_M = Matrix{Float64}(I, M, M)

        if P == 1
            A[:, 1:M] = Pi_mat + I_M
        else
            gamma1 = gamma_cols_for_lag(1)
            if any(c -> !(c in cn), gamma1)
                missing = filter(c -> !(c in cn), gamma1)
                throw(ArgumentError("Cannot map to VAR: missing lag-1 coefficient columns in VECM output: $(missing)."))
            end
            A[:, 1:M] = Pi_mat + I_M + Matrix(B_named[:, gamma1])

            for lag_idx in 2:(P - 1)
                g_curr = gamma_cols_for_lag(lag_idx)
                g_prev = gamma_cols_for_lag(lag_idx - 1)
                if any(c -> !(c in cn), g_curr)
                    missing = filter(c -> !(c in cn), g_curr)
                    throw(ArgumentError("Cannot map to VAR: missing lag-$(lag_idx) coefficient columns in VECM output: $(missing)."))
                end
                A[:, (M * (lag_idx - 1) + 1):(M * lag_idx)] = Matrix(B_named[:, g_curr]) - Matrix(B_named[:, g_prev])
            end

            g_last = gamma_cols_for_lag(P - 1)
            if any(c -> !(c in cn), g_last)
                missing = filter(c -> !(c in cn), g_last)
                throw(ArgumentError("Cannot map to VAR: missing lag-$(P - 1) coefficient columns in VECM output: $(missing)."))
            end
            A[:, (M * (P - 1) + 1):(M * P)] = -Matrix(B_named[:, g_last])
        end

        # Copy deterministic and exogenous coefficients unchanged.
        for c in det_exo_cols
            A[:, c] = B_named[:, c]
        end

        ret[:A_VAR] = A
    end

    return ret
end
