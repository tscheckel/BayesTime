# function to create companion form of VAR coefficient matrice
# author: tobias scheckel

# INPUTS:
# - A:              MxK-matrix, VAR coefficients (ATTENTION: coefficients associated with deterministic terms can be included but need to be ordered last)
# - M:              scalar, number of equations
# - P:              scalar, number of lags

# OUTPUTS:
# dictionary containing the following:
# - Acm:            (MP)x(MP)-matrix, VAR coefficients in companion form
# - J:              Mx(MP)-matrix

# ----- FUNCTION BODY -----
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