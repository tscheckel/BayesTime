"""
    VAR_OLS(data, Exo, P = 1, cons = true, quant = [.025; .5; .975])

Estimate a Vector Autoregression via Ordinary Least Squares. 

# INPUTS:
- `data::Matrix{Float64}`: Data matrix where rows are observations and columns are variables
- `Exo`: Optional matrix of exogenous variables (default: nothing)
- `P::Int`: Number of lags in the VAR (default: 1)
- `cons::Bool`: Whether to include a constant term (default: true)
- `quant::Vector{Float64}`: Quantiles for confidence intervals (default: [.

# OUTPUTS:
- `Dict` containing:
    - `:Alpha`: 3D array of VAR coefficients with dimensions (K, M, length(quant))
    - `:Sigma_hat`: Estimated covariance matrix of residuals
    - `:Y`: Dependent variable matrix
    - `:X`: Regressor matrix

# Author
Tobias Scheckel
"""
function VAR_OLS(;
    data::Matrix{Float64}, Exo = nothing,
    P = 1, cons = true,
    quant::Vector{Float64} = [.025; .5; .975]
)
    # get quantiles of standard Gaussian (for CIs)
    normq = quantile.(Normal(), quant)
    
    # ---- DATA ----
    Y = data[P+1:end,:]
    T, M = size(Y)
    X = mlag(data, P)
    if !isnothing(Exo)
        X = hcat(X, Exo)
    end
    if cons 
        X = hcat(X, fill(1.0, T))
    end
    K = size(X, 2)
    
    # ESTIMATION
    y = vec(Y)
    X_kron = kron(I(M), X)

    # OLS estimate of VAR coeffs:
    Alpha_hat = reshape(((X_kron'X_kron)\I)*(X_kron'y), K, M)
    # compute residuals
    U = Y - X*Alpha_hat
    # OLS estimate of 
    Sigma_hat = U'U/(T-K)
    
    # variance covariance matrix of Alpha_hat
    Alpha_V = kron(Sigma_hat, X'X)

    # for each confidence level: compute CIs of parameter and IRF estimates
    Alpha = fill(NaN, K, M, length(quant))
    for qx in 1:length(quant)
        Alpha[:,:,qx] = Alpha_hat + reshape(diag(Alpha_V), K, M)*normq[qx]
    end

    retDict = Dict(
        :Alpha => Alpha,
        :Sigma_hat => Sigma_hat,
        :Y => Y,
        :X => X
    )

    return retDict
end