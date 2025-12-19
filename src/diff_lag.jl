"""
    diff_lag(x, tlag)

Compute lagged differences of a time series vector.

# Arguments
- `x::AbstractVector`: Input time series vector
- `tlag::Int`: Number of lags for the difference operation

# Returns
Vector of length `length(x) - tlag` containing the lagged differences: 
`x[t] - x[t-tlag]` for t = tlag+1, ..., T.

# Example
```julia
x = [1.0, 2.0, 4.0, 7.0, 11.0]
diff_lag(x, 1)  # Returns [1.0, 2.0, 3.0, 4.0]
diff_lag(x, 2)  # Returns [3.0, 5.0, 7.0]
```

# Author
Tobias Scheckel
"""
function diff_lag(x::AbstractVector, tlag::Int)
    xdiff = x[tlag+1:end] - x[1:end-tlag]
    return xdiff
end