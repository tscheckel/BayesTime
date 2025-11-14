
x = 1e-8

function offsetting(x;
    thrshld1 = 1e-6, thrshld2 = 1e6)
    if x < thrshld1
        y = thrshld1
    elseif x > thrshld2
        y = thrshld2
    else
        y = copy(x)
    end
    
    return y
end

