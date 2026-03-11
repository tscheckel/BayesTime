"""
    get_companion(; A, M, P)

Create companion-form matrices for a VAR(P).

# Inputs
- `A::Matrix{Float64}`: `M x K` VAR coefficient matrix. Lag blocks must appear first as
  `[A_1 A_2 ... A_P]`. Deterministic/exogenous coefficients may be appended after lag blocks.
- `M::Int`: number of equations (endogenous variables)
- `P::Int`: number of lags

# Returns
- `Dict` with:
  - `"Acm"`: `(M*P) x (M*P)` companion matrix
  - `"J"`: `M x (M*P)` selector matrix

# Author
Tobias Scheckel
"""
function get_companion(;A::Matrix{Float64}, M::Int, P::Int)
    ## ----- Generate J -----
    J = hcat(I(M), zeros(M, M*(P-1)))
    # J[1:M,1:M] .= I(M)  # Set the first M x M block to an identity matrix

    ## ----- Generate -----
    Acm = zeros(M*P,M*P) # initialize
    # Fill in values
    for ll in 1:P
        # First line is coefficients for lags
        Acm[1:M,((ll-1)*M+1):(ll*M)] .= A[:,((ll-1)*M+1):(ll*M)]
        if ll > 1
            Acm[((ll-1)*M+1):(ll*M),((ll-2)*M+1):((ll-1)*M)] .= I(M)
        end
    end

    return Dict("Acm" => Acm, "J" => J)
end
