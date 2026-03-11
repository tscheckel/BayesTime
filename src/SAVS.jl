"""
    SAVS(beta, X)

Apply Signal-Adaptive Variable Selection (SAVS) shrinkage to coefficient estimates.

# Inputs
- `beta`: `K`-vector of parameter estimates
- `X`: `n x K` design matrix

# Returns
- `gamma`: `K`-vector of sparsified SAVS estimates

# Author
Tobias Scheckel
"""
function SAVS(beta, X)
    K = length(beta) # number of parameters
    gamma = zeros(K) # initialize output vector
    # sparsify each element in turn
    for j in 1:K
        κ = 1/beta[j]^2 # compute kappa
        Xj_sqnorm = sum(X[:, j].^2) # compute norm of design matrix column

        gamma[j] = sign(beta[j])*maximum([abs(beta[j])*Xj_sqnorm - κ, 0])/Xj_sqnorm # following Hauzenberger et al. (2021, JBES)
    end
    # return sparsified estimates
    return gamma
end
