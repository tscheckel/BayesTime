"""
    get_IRF(; A, shock, P, NHOR)

Compute impulse response functions for a VAR(P) to a given shock vector.

# Inputs
- `A::Matrix{Float64}`: `M x K` matrix of VAR coefficients
- `shock::Vector{Float64}`: `M`-vector with initial shock at horizon 0
- `P::Int`: number of VAR lags
- `NHOR::Int`: number of response horizons to compute

# Returns
- `IRF`: `(NHOR+1) x M` matrix of responses, where row 1 is the impact response (`h=0`)

# Author
Tobias Scheckel
"""
function get_IRF(;
    A::Matrix{Float64}, shock::Vector{Float64}, P::Int, NHOR::Int
    )
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
