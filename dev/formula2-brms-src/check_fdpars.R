function (x) 
{
    stopifnot(is.null(x) || is.list(x))
    pos_pars <- c("sigma", "shape", "nu", "phi", "kappa", "beta", 
        "disc", "bs", "ndt", "theta")
    prob_pars <- c("zi", "hu", "bias", "quantile")
    for (dp in names(x)) {
        apc <- dpar_class(dp)
        value <- x[[dp]]$value
        if (apc %in% pos_pars && value < 0) {
            stop2("Parameter '", dp, "' must be positive.")
        }
        if (apc %in% prob_pars && (value < 0 || value > 1)) {
            stop2("Parameter '", dp, "' must be between 0 and 1.")
        }
    }
    invisible(TRUE)
}
