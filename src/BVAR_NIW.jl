"""
    BVAR_NIW(; data, Exo = nothing, P = 1, cons = true, trend = false,
               lambda = .2, own_mean = 0.0, det_var = 100.0,
               Alpha_M0 = nothing, Alpha_V0 = nothing, Sigma_S0 = nothing, Sigma_nu0 = nothing,
               nsave = 10_000, stable = false, quant = [.025; .5; .975])

Estimate a Bayesian Vector Autoregression with a conjugate Normal-inverse-Wishart prior.

The model is Y = X*Alpha + U, rows of U iid N(0, Sigma), where X = [Y_{t-1} ... Y_{t-P} Exo trend const]
uses the same design matrix as `VAR_OLS` (built by `get_VAR_design`). The prior is

    Sigma ~ IW(Sigma_S0, Sigma_nu0)
    vec(Alpha) | Sigma ~ N(vec(Alpha_M0), Sigma ⊗ Alpha_V0)

so the posterior is again NIW with

    Alpha_V   = (Alpha_V0^-1 + X'X)^-1
    Alpha_M   = Alpha_V*(Alpha_V0^-1*Alpha_M0 + X'Y)
    Sigma_S   = Sigma_S0 + (Y - X*Alpha_M)'(Y - X*Alpha_M) + (Alpha_M - Alpha_M0)'Alpha_V0^-1(Alpha_M - Alpha_M0)
    Sigma_nu  = Sigma_nu0 + T

and draws are iid (no MCMC). By default, a Minnesota-type prior is used: prior mean `own_mean`
on the first own lag and zero elsewhere, prior variance (lambda/l)^2/sig2_j for lag l of
variable j and `det_var` for exogenous/deterministic terms, Sigma_S0 = diag(sig2)*(Sigma_nu0 - M - 1)
and Sigma_nu0 = M + 2, where sig2_j is the residual variance of a univariate AR(P) for variable j.
Because Sigma also scales the prior on Alpha (conjugacy), the Minnesota cross-variable scaling
sig2_i/sig2_j is implied. Any of `Alpha_M0`, `Alpha_V0` (KxK), `Sigma_S0`, `Sigma_nu0` can be
supplied to override the default.

# Arguments
- `data`: TxM data matrix or `DataFrame` (rows are observations)
- `Exo`: optional TxN matrix or `DataFrame` of exogenous variables (default: nothing)
- `P::Int`: number of lags (default: 1)
- `cons::Bool`: include a constant (default: true)
- `trend::Bool`: include a linear time trend (default: false)
- `lambda::Float64`: overall Minnesota tightness (default: .2)
- `own_mean::Float64`: prior mean of the first own lag; 1.0 for a random-walk prior (default: 0.0)
- `det_var::Float64`: prior variance of exogenous/deterministic coefficients (default: 100.0)
- `Alpha_M0`, `Alpha_V0`, `Sigma_S0`, `Sigma_nu0`: optional user-supplied prior moments
- `nsave::Int`: number of posterior draws (default: 10_000)
- `stable::Bool`: keep only draws with a stable companion matrix (default: false)
- `quant::Vector{Float64}`: posterior quantiles to report (default: [.025; .5; .975])

# Returns
- `Dict` containing:
    - `:Alpha`: KxMxlength(quant) posterior quantiles of the VAR coefficients
    - `:Sigma`: MxM posterior mean of the error covariance
    - `:Alpha_store`, `:Sigma_store`: KxMxnsave and MxMxnsave posterior draws
    - `:Alpha_M`, `:Alpha_V`, `:Sigma_S`, `:Sigma_nu`: posterior moments
    - `:Alpha_M0`, `:Alpha_V0`, `:Sigma_S0`, `:Sigma_nu0`: prior moments
    - `:Y`, `:X`: dependent variables and design matrix

# Author
Tobias Scheckel
"""
function BVAR_NIW(;
    data::Union{AbstractMatrix{<:Real}, DataFrame}, Exo::Union{AbstractMatrix{<:Real}, DataFrame, Nothing} = nothing,
    P::Int = 1, cons::Bool = true, trend::Bool = false,
    lambda::Float64 = .2, own_mean::Float64 = 0.0, det_var::Float64 = 100.0,
    Alpha_M0::Union{AbstractMatrix{<:Real}, Nothing} = nothing, Alpha_V0::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
    Sigma_S0::Union{AbstractMatrix{<:Real}, Nothing} = nothing, Sigma_nu0::Union{Real, Nothing} = nothing,
    nsave::Int = 10_000, stable::Bool = false,
    quant::Vector{Float64} = [.025; .5; .975]
)
    # ---- FOR DEBUGGING PURPOSES ONLY ----
    # data = randn(200, 3); Exo = nothing; P = 2; cons = true; trend = false
    # lambda = .2; own_mean = 0.0; det_var = 100.0
    # Alpha_M0 = nothing; Alpha_V0 = nothing; Sigma_S0 = nothing; Sigma_nu0 = nothing
    # nsave = 1_000; stable = false; quant = [.025; .5; .975]

    # ---- DATA ----
    # same design matrix as VAR_OLS, so both estimators line up row for row
    design = get_VAR_design(data = data, Exo = Exo, P = P, cons = cons, trend = trend)
    Y = design[:Y]
    X = design[:X]
    reg_names = design[:reg_names]
    data_var_names = design[:var_names]
    T, M = size(Y)
    K = size(X, 2)
    Klag = M*P

    # ---- PRIOR ----
    # univariate AR(P) residual variances (same estimation sample as the VAR) for Minnesota scaling;
    # the own lags of variable j sit in columns j, M+j, ..., (P-1)*M+j of X
    sig2 = fill(NaN, M)
    for j in 1:M
        X_j = hcat(X[:,j:M:Klag], fill(1.0, T))
        u_j = Y[:,j] - X_j*(X_j\Y[:,j])
        sig2[j] = u_j'u_j/(T-P-1)
    end

    if isnothing(Alpha_M0)
        Alpha_M0 = zeros(K, M)
        Alpha_M0[1:M,1:M] = own_mean*I(M)
    end
    if isnothing(Alpha_V0)
        # conjugacy forces a Kronecker structure, so the prior variance can depend on the
        # regressor (lag l, variable j) but not on the equation
        Alpha_V0_diag = fill(det_var, K)
        for l in 1:P, j in 1:M
            Alpha_V0_diag[(l-1)*M+j] = (lambda/l)^2/sig2[j]
        end
        Alpha_V0 = Diagonal(Alpha_V0_diag)
    end
    if isnothing(Sigma_nu0)
        Sigma_nu0 = M + 2
    end
    if isnothing(Sigma_S0)
        # chosen so that the prior mean of Sigma equals diag(sig2)
        Sigma_S0 = Diagonal(sig2)*(Sigma_nu0 - M - 1)
    end
    Alpha_M0 = Matrix{Float64}(Alpha_M0)
    Alpha_V0 = Matrix{Float64}(Alpha_V0)
    Sigma_S0 = Matrix{Float64}(Sigma_S0)

    # ---- POSTERIOR ----
    Alpha_V0_inv = Alpha_V0\I
    Alpha_V = forcePD(Matrix(Symmetric((Alpha_V0_inv + X'X)\I)))
    Alpha_M = Alpha_V*(Alpha_V0_inv*Alpha_M0 + X'Y)
    U = Y - X*Alpha_M
    Sigma_S = forcePD(Sigma_S0 + U'U + (Alpha_M - Alpha_M0)'Alpha_V0_inv*(Alpha_M - Alpha_M0))
    Sigma_nu = Sigma_nu0 + T

    # ---- POSTERIOR DRAWS ----
    Alpha_store = fill(NaN, K, M, nsave)
    Sigma_store = fill(NaN, M, M, nsave)
    Alpha_V_chol = cholesky(Symmetric(Alpha_V)).L
    Sigma_post = InverseWishart(Sigma_nu, Sigma_S)

    # with stable = true, unstable draws are rejected; cap the attempts so a posterior that sits
    # (almost) entirely outside the stationary region does not loop forever
    max_tries = stable ? 100*nsave : nsave
    save_id = 0
    ntries = 0
    while (save_id < nsave) & (ntries < max_tries)
        ntries += 1
        Sigma_draw = rand(Sigma_post)
        # vec(Alpha) ~ N(vec(Alpha_M), Sigma ⊗ Alpha_V)  <=>  Alpha = Alpha_M + L_V*Z*L_Sigma'
        Alpha_draw = Alpha_M + Alpha_V_chol*randn(K, M)*cholesky(Symmetric(Sigma_draw)).U

        # CHECK STABILITY
        if stable
            Acm = get_companion(A = Matrix(Alpha_draw[1:Klag,:]'), M = M, P = P)["Acm"]
            if maximum(abs.(eigvals(Acm))) >= 1.0
                continue
            end
        end

        # STORE DRAWS
        save_id += 1
        Alpha_store[:,:,save_id] = Alpha_draw
        Sigma_store[:,:,save_id] = Sigma_draw
    end
    if save_id < nsave
        @warn "Only $save_id of $nsave stable draws found in $max_tries attempts."
        Alpha_store = Alpha_store[:,:,1:save_id]
        Sigma_store = Sigma_store[:,:,1:save_id]
    elseif stable
        println("Stability acceptance rate: $(round(nsave/ntries, digits = 3))")
    end

    # ---- POSTERIOR SUMMARIES ----
    Alpha = fill(NaN, K, M, length(quant))
    for i in 1:K, j in 1:M
        Alpha[i,j,:] = quantile(Alpha_store[i,j,:], quant)
    end
    # without stability truncation the posterior mean is analytic; otherwise use the draws
    Sigma = stable ? dropdims(sum(Sigma_store, dims = 3), dims = 3)/size(Sigma_store, 3) : Sigma_S/(Sigma_nu - M - 1)

    Alpha = NamedArray(
        Alpha,
        (reg_names, data_var_names, quant),
        ("regressors", "vars", "quantiles")
    )
    Sigma = NamedArray(
        Sigma,
        (data_var_names, data_var_names),
        ("vars", "vars")
    )

    retDict = Dict(
        # posterior summaries
        :Alpha => Alpha,
        :Sigma => Sigma,
        # posterior draws
        :Alpha_store => Alpha_store,
        :Sigma_store => Sigma_store,
        # posterior moments
        :Alpha_M => Alpha_M,
        :Alpha_V => Alpha_V,
        :Sigma_S => Sigma_S,
        :Sigma_nu => Sigma_nu,
        # prior moments
        :Alpha_M0 => Alpha_M0,
        :Alpha_V0 => Alpha_V0,
        :Sigma_S0 => Sigma_S0,
        :Sigma_nu0 => Sigma_nu0,
        # data
        :Y => Y,
        :X => X
    )

    return retDict
end
