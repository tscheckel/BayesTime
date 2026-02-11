"""
    VAR_OLS(data, Exo, P = 1, cons = true, quant = [.025; .5; .975])

Estimate a Vector Autoregression via Ordinary Least Squares. 

# INPUTS:
- `data::Matrix{Float64}`: TxM-dimensional Data matrix where rows are observations and columns are variables
- `Exo`: TxM-dimensional Optional matrix of exogenous variables (default: nothing)
- `P::Int`: Number of lags in the VAR (default: 1)
- `cons::Bool`: Whether to include a constant term (default: true)
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
    P::Int = 1, cons::Bool = true,
    quant::Vector{Float64} = [.025; .5; .975]
)
    # CHECK INPUT DATATYPES
    if isa(data, DataFrame)
        data = NamedArray(
            Matrix{Float64}(data),
            (collect(axes(data, 1)), collect(names(data))),
            ("obs", "vars")
        )
    end 
    if isa(Exo, DataFrame)
        Exo = NamedArray(
            Matrix{Float64}(Exo),
            (collect(axes(Exo, 1)), collect(names(Exo))),
            ("obs", "exo")
        )
    end
    # get quantiles of standard Gaussian (for CIs)
    normq = quantile.(Normal(), quant)

    data_var_names = collect(1:size(data, 2))
    if data isa NamedArray
        data_var_names = collect(names(data, 2))
    end
    
    # ---- DATA ----
    # create matrix of dependent variables and regressors
    Y = Matrix{Float64}(data[P+1:end,:])
    T, M = size(Y)
    # create design matrix
    X_lag = mlag(data, P)
    X = Matrix{Float64}(X_lag)
    reg_names = collect(1:size(X, 2))
    if X_lag isa NamedArray
        reg_names = collect(names(X_lag, 2))
    end
    # add exogenous regressors if provided
    if !isnothing(Exo)
        Exo_slice = Exo[P+1:end,:]
        X = hcat(X, Matrix{Float64}(Exo_slice))
        exo_names = ["exo_$i" for i in 1:size(Exo_slice, 2)]
        if Exo_slice isa NamedArray
            exo_names = collect(names(Exo_slice, 2))
        end
        reg_names = vcat(reg_names, exo_names)
    end
    if cons 
        X = hcat(X, fill(1.0, T))
        reg_names = vcat(reg_names, ["const"])
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
