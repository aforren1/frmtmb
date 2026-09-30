plot.brmsfit <- 
function (x, pars = NA, combo = c("hist", "trace"), nvariables = 5, 
    N = NULL, variable = NULL, regex = FALSE, fixed = FALSE, 
    bins = 30, theme = NULL, plot = TRUE, ask = TRUE, newpage = TRUE, 
    ...) 
{
    contains_draws(x)
    nvariables <- use_alias(nvariables, N)
    if (!is_wholenumber(nvariables) || nvariables < 1) {
        stop2("Argument 'nvariables' must be a positive integer.")
    }
    variable <- use_variable_alias(variable, x, pars, fixed = fixed)
    if (is.null(variable)) {
        variable <- default_plot_variables(x)
        regex <- TRUE
    }
    draws <- as.array(x, variable = variable, regex = regex)
    variables <- dimnames(draws)[[3]]
    if (!length(variables)) {
        stop2("No valid variables selected.")
    }
    if (plot) {
        default_ask <- devAskNewPage()
        on.exit(devAskNewPage(default_ask))
        devAskNewPage(ask = FALSE)
    }
    n_plots <- ceiling(length(variables)/nvariables)
    plots <- vector(mode = "list", length = n_plots)
    for (i in seq_len(n_plots)) {
        sub <- ((i - 1) * nvariables + 1):min(i * nvariables, 
            length(variables))
        sub_vars <- variables[sub]
        sub_draws <- draws[, , sub_vars, drop = FALSE]
        plots[[i]] <- bayesplot::mcmc_combo(sub_draws, combo = combo, 
            bins = bins, gg_theme = theme, ...)
        if (plot) {
            plot(plots[[i]], newpage = newpage || i > 1)
            if (i == 1) {
                devAskNewPage(ask = ask)
            }
        }
    }
    invisible(plots)
}
