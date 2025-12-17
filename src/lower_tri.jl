"""
    lower_tri(X::Matrix; diag::Bool = false)

Create vector of lower triangular elements of a square matrix.

# Arguments
- `X::Matrix`: Square matrix from which to extract lower triangular elements
- `diag::Bool`: Whether to include the diagonal elements (default: false)

# Returns
Vector containing the lower triangular elements of the matrix.

# Author
Tobias Scheckel
"""
lower_tri = function(X::Matrix; diag::Bool = false)
    if diag
        lt = X[tril(trues(size(X)))]
    else
        lt = X[tril(trues(size(X)), -1)]
    end
    return lt
end