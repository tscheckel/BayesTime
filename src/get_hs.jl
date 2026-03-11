"""
    get_hs(; bdraw, λ, τ, ν, ζ)

Draw/update Horseshoe prior parameters.

# Inputs
- `bdraw::Vector{Float64}`: regression parameter draw
- `λ::Vector{Float64}`: local variance components
- `τ::Float64`: global variance component
- `ν::Vector{Float64}`: local auxiliary variables
- `ζ::Float64`: global auxiliary variable

# Returns
- `Dict` with:
  - `"ψ"`: prior variances (`λ .* τ`)
  - `"λ"`: updated local variance components
  - `"τ"`: updated global variance component
  - `"ν"`: updated local auxiliary variables
  - `"ζ"`: updated global auxiliary variable

# Author
Tobias Scheckel
"""
function get_hs(;
    bdraw::Vector{Float64},
    λ::Vector{Float64}, τ::Float64,
    ν::Vector{Float64}, ζ::Float64
    )
    
    k = length(bdraw)
    
    # global & local prior variance components
    for j in 1:k
        λ[j] = rand(InverseGamma(1.0, 1.0/ν[j] + bdraw[j]^2/(2.0 * τ)))
    end
    τ = rand(InverseGamma((k+1.0)/2.0, 1.0/ζ + sum(bdraw.^2 ./ λ)/2.0))
    
    # auxiliary variables
    for j in 1:k
        ν[j] = rand(InverseGamma(1.0, 1.0 + 1.0/λ[j]))
    end
    ζ = rand(InverseGamma(1.0, 1.0 + 1.0 / τ))
    
    ret = Dict(
        "ψ" => (λ * τ),
        "λ" => λ,
        "τ" => τ,
        "ν" => ν,
        "ζ" => ζ
    )
    return ret
end
