# function to evaluate convergence of a MCMC-series

# author: tobias scheckel

# inputs:
# - series:   vector, MCMC series
# - nburn:    scalar, number of burn-in draws
# - title:    string, plot title

using MCMCDiagnosticTools
using Plots

function eval_conv(series, nburn, titile = "")
    p_geweke = round(gewekediag(series)[2],digits = 4)
    plot(series, title = "p-value Geweke test: $p_geweke",
    color = :grey)
    vline!([nburn], linestyle = :dash, linecolor = :red)
end
