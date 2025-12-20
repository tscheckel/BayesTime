"""
    data_trans(; data, tcodes, tlag = 1)

Apply data transformations to time series variables according to transformation codes.

# Arguments
- `data::Union{DataFrame, AbstractMatrix}`: TxM matrix where rows are time periods and columns are variables
- `tcodes::Vector{Int}`: M-vector of transformation codes for each variable
- `tlag::Int`: Number of periods for lag/difference operations (default: 1)

# Transformation Codes
- `1`: Level (no transformation)
- `2`: First difference
- `3`: Second difference (not yet implemented)
- `4`: Natural logarithm
- `5`: First difference of natural logarithm (× 100)
- `6`: Second difference of natural logarithm (not yet implemented)

# Returns
TxM matrix with transformed variables. Note that differencing operations create NaN values
at the beginning of the series.

# Author
Tobias Scheckel
"""
function data_trans(;
    data::Union{DataFrame, AbstractMatrix}, tcodes::Vector{Int}, tlag::Int = 1)
    
    # get dimensions
    T, M = size(data)
    
    # initialize output object
    data_out = copy(data) #fill(NaN, T, M)
    data_out .= NaN # set all values to NaN (as differencing will create NaNs at start)
    
    # loop over variables and apply transformations
    for i in 1:M
        if tcodes[i] == 1
            # level data, no transformation
            data_out[:, i] = data[:, i]
        elseif tcodes[i] == 2
            # first difference
            data_out[tlag+1:end, i] = diff_lag(data[:, i], tlag)
        elseif tcodes[i] == 3
            # second difference
            data_out[tlag*2+1:end, i] = diff_lag(diff_lag(data[:, i], tlag), tlag)
        elseif tcodes[i] == 4
            # log level
            data_out[:, i] = log.(data[:, i])
        elseif tcodes[i] == 5
            # log first difference
            data_out[tlag+1:end, i] = diff_lag(log.(data[:,i]), tlag)*100
        elseif tcodes[i] == 6
            # log second difference
            data_out[tlag*2+1:end, i] = diff_lag(diff_lag(log.(data[:, i]), tlag), tlag)*100
        else
            error("Unknown transformation code: $(tcodes[i])")
        end
    end
    
    return data_out
end