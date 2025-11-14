# function to implement the efficient estimation algorithm for large VARs with non-conjugate priors and SV
# author: tobias scheckel
# created: May 16 2025

# INPUTS:
# Y:        TxM-matrix of dependent variables
# X:        TxK-matrix of regressors
# Π:        KxM-matrix of VAR coefficients (previous draw)
# A:        MxM-matrix governing contemporaneous interactions (should be upper-triangular)
# Λ:        TxM-matrix of VAR residual variances
# Π_V0:     KxM-matrix of VAR coefficient prior variances

# OUTPUTS:
# Π:        KxM-matrix of VAR coefficients (updated draw)

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