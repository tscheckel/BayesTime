function offsetting(x;
    thrshld1 = 1e-6, thrshld2 = 1e6)
    y = copy(x)
    y[x < thrshld1] = thrshld1
    y[x > thrshld2] = thrshld2
    return y
end