"""
    VECM(; data, Exo=nothing, P::Int=1, r::Int, cons::Bool=true, trend::Bool=false,
          estim::String="ML", map_to_var::Bool=false)

Wrapper around `tsDyn::VECM` via `RCall`.

# Inputs
- `data`: T x M matrix-like object (`AbstractMatrix` or `DataFrame`)
- `Exo`: optional T x K exogenous matrix-like object (`nothing` by default)
- `P`: VAR lag order in levels. Internally mapped to `tsDyn::VECM(lag = P - 1)`.
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
    
    # Capture variable names before matrix conversion (plain matrices do not carry names).
    data_col_names = if data isa DataFrame
        String.(names(data))
    elseif data isa NamedArray
        String.(names(data, 2))
    else
        nothing
    end

    exo_col_names = if Exo isa DataFrame
        String.(names(Exo))
    elseif Exo isa NamedArray
        String.(names(Exo, 2))
    else
        nothing
    end

    # Normalize accepted input types (DataFrame/Matrix-like) to dense Float64 matrices.
    # This guarantees predictable element types before passing data to R.
    data_mat = Matrix{Float64}(data)
    T, M = size(data_mat)

    if !isnothing(Exo)
        exo_mat = Matrix{Float64}(Exo)
        size(exo_mat, 1) == T || throw(ArgumentError("Exo must have the same number of rows as data."))
    else
        exo_mat = nothing
    end

    # Match tsDyn::VECM deterministic-term options.
    include = if trend && cons
        "both"
    elseif trend
        "trend"
    elseif cons
        "const"
    else
        "none"
    end

    # Push Julia values into the R session; names are reused in the R block below.
    @rput data_mat
    @rput exo_mat
    @rput P
    @rput r
    @rput include
    @rput estim
    @rput data_col_names
    @rput exo_col_names

    R"""
    library(tsDyn)

    if (!is.null(data_col_names)) {
      colnames(data_mat) <- data_col_names
    }
    if (!is.null(exo_mat) && !is.null(exo_col_names)) {
      colnames(exo_mat) <- exo_col_names
    }

    if (is.null(exo_mat)) {
      vecm <- VECM(data_mat, lag = P-1, r = r, include = include, estim = estim)
    } else {
      vecm <- VECM(data_mat, lag = P-1, r = r, include = include, estim = estim, exogen = exo_mat)
    }

    B <- coef(vecm)
    Pi <- coefPI(vecm)
    cn <- colnames(B)
    rn <- rownames(B)
    """

    # Pull estimated objects back from R into Julia.
    @rget B
    @rget Pi
    @rget cn
    @rget rn

    # R may return string-like vectors; normalize to Vector{String}.
    cn = String.(cn)
    rn = String.(rn)

    B_mat = Matrix{Float64}(B)
    Pi_mat = Matrix{Float64}(Pi)

    # Attach row/column labels so downstream indexing stays name-based.
    B_named = NamedArray(B_mat, (rn, cn), ("Eq", "Coeff"))

    # Prefer user-provided data names; otherwise fall back to names from Pi/generic labels.
    var_names = if !isnothing(data_col_names) && length(data_col_names) == M
        data_col_names
    elseif size(Pi_mat, 2) == M && Pi isa NamedArray
        String.(names(Pi, 2))
    else
        ["Var" * string(i) for i in 1:M]
    end

    # Keep equation names from R rows, but enforce economically meaningful variable names on columns.
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
        # tsDyn names short-run lag-difference terms as "<var> -<lag>".
        function gamma_cols_for_lag(lag_idx::Int)
            return [string(v, " -", lag_idx) for v in var_names]
        end

        # Deterministic/exogenous terms that are not ECT and not lagged endogenous terms.
        lag_pattern = r" -[0-9]+$"
        det_exo_cols = [c for c in cn if c != "ECT" && !occursin(lag_pattern, c)]

        # Target VAR coefficient layout:
        # [A1 A2 ... AP | deterministic/exogenous terms]
        A_col_lag = [string(v, "_L", lag_idx) for lag_idx in 1:P for v in var_names]
        A_colnames = vcat(A_col_lag, det_exo_cols)
        A = NamedArray(
            fill(NaN, M, length(A_colnames)),
            (rn, A_colnames),
            ("Eq", "Coeff")
        )

        I_M = Matrix{Float64}(I, M, M)

        if P == 1
            # For VECM with one lag:
            # Δy_t = Π y_{t-1} + ...
            # implies y_t = (I + Π) y_{t-1} + ...
            A[:, 1:M] = Pi_mat + I_M
        else
            gamma1 = gamma_cols_for_lag(1)
            if any(c -> !(c in cn), gamma1)
                missing = filter(c -> !(c in cn), gamma1)
                throw(ArgumentError("Cannot map to VAR: missing lag-1 coefficient columns in VECM output: $(missing)."))
            end
            # A1 = I + Π + Γ1
            A[:, 1:M] = Pi_mat + I_M + Matrix(B_named[:, gamma1])

            for lag_idx in 2:(P - 1)
                g_curr = gamma_cols_for_lag(lag_idx)
                g_prev = gamma_cols_for_lag(lag_idx - 1)
                if any(c -> !(c in cn), g_curr)
                    missing = filter(c -> !(c in cn), g_curr)
                    throw(ArgumentError("Cannot map to VAR: missing lag-$(lag_idx) coefficient columns in VECM output: $(missing)."))
                end
                # Aj = Γj - Γ{j-1}, for j = 2, ..., P-1
                A[:, (M * (lag_idx - 1) + 1):(M * lag_idx)] = Matrix(B_named[:, g_curr]) - Matrix(B_named[:, g_prev])
            end

            g_last = gamma_cols_for_lag(P - 1)
            if any(c -> !(c in cn), g_last)
                missing = filter(c -> !(c in cn), g_last)
                throw(ArgumentError("Cannot map to VAR: missing lag-$(P - 1) coefficient columns in VECM output: $(missing)."))
            end
            # AP = -Γ{P-1}
            A[:, (M * (P - 1) + 1):(M * P)] = -Matrix(B_named[:, g_last])
        end

        # Deterministic/exogenous coefficients are already in levels form, so copy directly.
        for c in det_exo_cols
            A[:, c] = B_named[:, c]
        end

        # Optional output; base return still includes B/Pi even when mapping is skipped.
        ret[:A_VAR] = A
    end

    return ret
end
