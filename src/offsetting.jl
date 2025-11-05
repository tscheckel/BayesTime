function offsetting(x;
    thrshld1 = 1e-6, thrshld2 = 1e6)
    # if x is float, convert to vector (otherwise subsetting does not work)
    if typeof(x) == Float64
        x = [x]
    end
    y = copy(x)
    y[x .< thrshld1] .= thrshld1
    y[x .> thrshld2] .= thrshld2
    return y
end