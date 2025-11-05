# function to create matrix of lagged observations
# author: tobias scheckel

# INPUTS:
# - X:              TxM-matrix (Traw = number of obs, M number of variables), data
# - P:              scalar, number of lags

# OUTPUT:
# - Xlag:           (T-P)x(M*P)-matrix, lagged variables

# ----- FUNCTION BODY -----
function mlag(X::Matrix{Float64}, P::Int)
    # X = convert(Matrix{Float64}, X)  # Ensure X is a Float64 matrix
    Traw, N = size(X)
    Xlag = fill(NaN, Traw-P, P * N)  # Initialize Xlag with NaNs
    
    for ii in 1:P
        Xlag[:, (N * (ii - 1) + 1):(N * ii)] = X[(P + 1 - ii):(Traw - ii), 1:N]
    end
    
    return Xlag
end