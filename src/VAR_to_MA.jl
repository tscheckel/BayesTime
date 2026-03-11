"""
    VAR_to_MA(A, M, P, Q)

Map VAR coefficients to moving-average coefficients.

# Inputs
- `A`: `M x (M*P)` matrix of VAR lag coefficients with blocks ordered as
  `[A_1 A_2 ... A_P]`
- `M`: number of endogenous variables
- `P`: number of VAR lags
- `Q`: number of MA lags to compute beyond `Psi_0`

# Returns
- `Psi`: `M x (M*(Q+1))` matrix containing MA blocks `[Psi_0 Psi_1 ... Psi_Q]`

# Notes
- Output includes `Psi_0 = I_M` as the first `M x M` block.
- `A` is assumed to contain only the stochastic VAR lag coefficients.
"""
function VAR_to_MA(
    A::AbstractMatrix{<:Real}, M::Int, P::Int, Q::Int
)
    # ----- Basic argument checks -----
    # Dimensions must define a valid VAR(P) system with at least one variable/lag/horizon.
    M < 1 && throw(ArgumentError("M must be >= 1."))
    P < 1 && throw(ArgumentError("P must be >= 1."))
    Q < 1 && throw(ArgumentError("Q must be >= 1."))

    # A must contain exactly the P lag blocks [A_1 ... A_P], each of size M x M.
    size(A, 1) != M && throw(ArgumentError("A must have M rows. Got size(A,1)=$(size(A,1)) and M=$(M)."))
    size(A, 2) != M*P && throw(ArgumentError("A must have M*P columns. Got size(A,2)=$(size(A,2)) and M*P=$(M*P)."))

    # Storage for [Psi_0, Psi_1, ..., Psi_Q] in block-column format.
    # Block q+1 (1-based) corresponds to Psi_q.
    Psi = fill(NaN, M, M*(Q+1))

    # Initialize zero-horizon MA coefficient: Psi_0 = I_M.
    Psi[:, 1:M] .= I(M)

    # Recursive MA mapping:
    # Psi_q = sum_{l=1}^{min(P,q)} A_l * Psi_{q-l}
    for q in 1:Q
        # Accumulator for current horizon q.
        Psi_q = zeros(M, M)
        for ll in 1:min(P, q)
            # Extract lag-l coefficient block A_l from A:
            # columns ((l-1)M+1) : (lM)
            A_ll = A[:, ((ll-1)*M+1):(ll*M)]

            # Extract Psi_{q-l} from block storage:
            # block index (q-l+1) in 1-based indexing.
            Psi_lag = Psi[:, ((q-ll)*M+1):((q-ll+1)*M)]

            # Add contribution of lag l to current MA block.
            Psi_q += A_ll*Psi_lag
        end

        # Write Psi_q into block (q+1), i.e., columns (qM+1):((q+1)M).
        Psi[:, (q*M+1):((q+1)*M)] .= Psi_q
    end

    return Psi
end
