# function to compute impulse response function of a VAR to a given shock
# author: tobias Scheckel

# INPUTS:
# - A:              MxK matrix of VAR coeffs
# - shock:          M-vector of shock
# - P:              scalar, number of VAR lags
# - NHOR:           scalar, number of impulse response horizons

# OUTPUTS:
# IRF:              (NHOR+1)xM-matrix of impulse reponse functions for each variable

# ---- FUNCTION BODY ----
function get_IRF(;A::Matrix{Float64}, shock::Vector{Float64}, P::Int, NHOR::Int)
    # get dimensions
    M = size(shock, 1)

    get_Cm = get_companion(A = A, M = M, P = P)
    J = get_Cm["J"] 
    Acm = get_Cm["Acm"]

    # storage
    IRF = fill(NaN, NHOR+1, M)
    IRF[1,:] = shock # initial response is equal to shock
    
    Cmi = Acm # intialize
    # compute response for each horizon
    for nhor in 1:NHOR
        IRF[nhor+1,:,] = J*Cmi*J'*shock
        Cmi = Cmi*Acm
    end
    
    return IRF
end