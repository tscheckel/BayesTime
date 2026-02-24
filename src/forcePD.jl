# FUNCTION TO FORCE A SQUARE MATRIX TO BE POSITIVE DEFINITE
# author: tobias scheckel

# INPUTS:
# X:            Matrix
# tol:          float, numerical tolerance

# OUTPUT:
# X_PD:         amended X, which is symmetric and has eigenvalues > 0

# ---- FUNCTION BODY ----
function forcePD(X::AbstractMatrix{<:Real}, tol::Float64=1e-12)
    # make symmertric
    if X != X'
        Xsimm = (X + X') ./ 2
    else
        Xsimm = copy(X)
    end

    # force all eigenvalues to be positive
    λ = eigen(Xsimm).values
    λmin = minimum(λ)
    δ = max(0.0, tol - λmin)
    X_PD = Xsimm + δ*I

    return X_PD
end