"""
    eval_conv(series, nburn, title="")

Produce trace plot and evaluate convergence of a MCMC series
usingthe Geweke diagnostic test.

# Arguments
- `series::Vector`: MCMC series
- `nburn::Number`: Number of burn-in draws
- `title::String`: Plot title (optional)

# Returns
Trace plot of the MCMC series with Geweke test p-value.

# Author
Tobias Scheckel
"""
function trace_plot(series::Vector, nburn::Int, title = "")
    p_geweke = round(gewekediag(series)[2],digits = 4)
    plot(series, title = "p-value Geweke test: $p_geweke",
    color = :grey)
    vline!([nburn], linestyle = :dash, linecolor = :red)
    display(current())
end
