"""
    mlag(X, P)

Create a matrix of lagged observations.

# Inputs
- `X::AbstractMatrix{<:Real}`: `T x M` data matrix
- `P::Int`: number of lags

# Returns
- `Xlag`: `(T-P) x (M*P)` matrix of lagged variables, ordered as
  `[X_{t-1} X_{t-2} ... X_{t-P}]`

# Notes
- If `X` is a `NamedArray`, names are propagated to the output.

# Author
Tobias Scheckel
"""
function mlag(X::AbstractMatrix{<:Real}, P::Int)
    X_raw = X
    X = Matrix{Float64}(X_raw)
    Traw, N = size(X)
    Xlag = fill(NaN, Traw-P, P * N)  # Initialize Xlag with NaNs
    
    for ii in 1:P
        Xlag[:, (N * (ii - 1) + 1):(N * ii)] = X[(P + 1 - ii):(Traw - ii), 1:N]
    end

    if X_raw isa NamedArray
        row_names_in = collect(names(X_raw, 1))
        col_names_in = collect(names(X_raw, 2))
        lag_col_names = [string(var_name, "_L", lag) for lag in 1:P for var_name in col_names_in]
        return NamedArray(
            Xlag,
            (row_names_in[P+1:end], lag_col_names),
            (dimnames(X_raw)[1], dimnames(X_raw)[2])
        )
    end

    return Xlag
end
