"""
    VAR_OLS(data, Exo, P = 1, cons = true, trend = false, quant = [.025; .5; .975])

Estimate a Vector Autoregression via Ordinary Least Squares. 

# INPUTS:
- `data::Matrix{Float64}`: TxM-dimensional Data matrix where rows are observations and columns are variables
- `Exo`: TxM-dimensional Optional matrix of exogenous variables (default: nothing)
- `P::Int`: Number of lags in the VAR (default: 1)
- `cons::Bool`: Whether to include a constant term (default: true)
- `trend::Bool`: Whether to include a deterministic linear time trend (default: false)
- `quant::Vector{Float64}`: Quantiles for confidence intervals (default: [.025; .5; .975])

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
    data::Union{AbstractMatrix{<:Real}, DataFrame}, Exo::Union{AbstractMatrix{<:Real}, DataFrame, Nothing} = nothing,
    P::Int = 1, cons::Bool = true, trend::Bool = false,
    quant::Vector{Float64} = [.025; .5; .975]
)
    # get quantiles of standard Gaussian (for CIs)
    normq = quantile.(Normal(), quant)

    # ---- DATA ----
    design = get_VAR_design(data = data, Exo = Exo, P = P, cons = cons, trend = trend)
    Y = design[:Y]
    X = design[:X]
    reg_names = design[:reg_names]
    data_var_names = design[:var_names]
    T, M = size(Y)
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
    
    # variance covariance matrix of vec(Alpha_hat)
    Alpha_V = kron(Sigma_hat, (X'X)\I)

    # for each confidence level: compute CIs of parameter and IRF estimates
    Alpha = fill(NaN, K, M, length(quant))
    for qx in 1:length(quant)
        Alpha[:,:,qx] = Alpha_hat + sqrt.(reshape(diag(Alpha_V), K, M))*normq[qx]
    end

    Alpha = NamedArray(
        Alpha,
        (reg_names, data_var_names, quant),
        ("regressors", "vars", "quantiles")
    )
    Sigma_hat = NamedArray(
        Sigma_hat,
        (data_var_names, data_var_names),
        ("vars", "vars")
    )

    retDict = Dict(
        :Alpha => Alpha,
        :Sigma_hat => Sigma_hat,
        :Y => Y,
        :X => X
    )

    return retDict
end
