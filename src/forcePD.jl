"""
    forcePD(X::AbstractMatrix{<:Real}, tol::Real=1e-12; maxiter::Int=12, growth::Real=10)

Return a symmetric positive-definite approximation of `X`.

For floating-point matrices, the function first attempts an eigenvalue-based
minimal diagonal shift. If that path fails, or for non-floating element types
(e.g. autodiff dual numbers), it falls back to a Cholesky-based strategy.

# Arguments
- `X`: Square real-valued matrix.
- `tol`: Minimum diagonal regularization.
- `maxiter`: Maximum number of Cholesky-jitter retries.
- `growth`: Multiplicative jitter growth factor per retry.

# Returns
- A symmetrized matrix adjusted to be positive definite when possible.

# Throws
- `ErrorException` if `X` contains `NaN` or `Inf`.
- Any unexpected factorization error is rethrown.
"""
function forcePD(X::AbstractMatrix{<:Real}, tol::Real=1e-12; maxiter::Int=12, growth::Real=10)
    if any(isnan.(X)) || any(isinf.(X))
        error("Input matrix contains NaN or Inf values")
    end

    # Work on an explicit symmetric copy to avoid mutating caller-owned storage.
    X_pd = Matrix((X + X') ./ 2)
    T = eltype(X_pd)

    # Prefer a one-shot minimal shift for standard floating-point matrices.
    if T <: AbstractFloat
        try
            λmin = minimum(eigvals(Symmetric(X_pd)))
            δ = max(zero(T), T(tol) - T(λmin))
            X_pd .+= δ .* I
            return X_pd
        catch
            # Fall through to Cholesky-jitter path.
        end
    end

    jitter = T(tol)
    for _ in 1:maxiter
        try
            cholesky(Symmetric(X_pd); check=true)
            return X_pd
        catch e
            if e isa PosDefException || e isa ZeroPivotException || e isa SingularException
                # Increase diagonal regularization geometrically until Cholesky succeeds.
                X_pd .+= jitter .* I
                jitter *= T(growth)
            else
                rethrow(e)
            end
        end
    end

    return X_pd
end
