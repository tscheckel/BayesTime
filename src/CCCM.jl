"""
    CCCM(; Y, X, Π, A, Λ, Π_M0, Π_V0)

Implement the equation-by-equation CCCM update for VAR coefficients under
stochastic volatility and non-conjugate priors.

# Inputs
- `Y`: `T x M` matrix of dependent variables
- `X`: `T x K` matrix of regressors
- `Π`: `K x M` matrix of VAR coefficients (previous draw)
- `A`: `M x M` matrix governing contemporaneous interactions (typically upper-triangular)
- `Λ`: `T x M` matrix of residual variances
- `Π_M0`: `K x M` matrix of prior means for VAR coefficients
- `Π_V0`: `K x M` matrix of prior variances for VAR coefficients

# Returns
- `Π`: updated `K x M` matrix of VAR coefficients

# Author
Tobias Scheckel
"""
function CCCM(;Y, X, Π, A, Λ, Π_M0, Π_V0)
    # get dimensions
    K, N = size(Π)
    # draw VAR coeffs eq-by-eq
    for i in 1:N
        Π_i = copy(Π)
        Π_i[:,i] .= 0
        rscl = vec(Λ[:,i:N].^0.5) # following Josh Chan's notation
        Y_i = vec((Y-X*Π_i)*A[i:N,:]')./rscl #eq 17
        X_i = kron(A[i:N,i],X)./rscl
        
        # post variance
        Π_i_V = Symmetric((Diagonal(1 ./Π_V0[:,i]) + X_i'X_i)\I)
        # post mean:
        if all(Π_M0[:,i] .== 0.0)
            Π_i_mean = Π_i_V*(X_i'Y_i)
        else
            Π_i_mean = Π_i_V*(X_i'Y_i + Diagonal(1 ./Π_V0[:,i])*Π_M0[:,i])
        end
        # sample ith column of A
        Π[:,i] = Π_i_mean + cholesky(Π_i_V).L*randn(K)
    end
    return Π
end
