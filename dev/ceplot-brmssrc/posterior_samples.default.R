posterior_samples.default <- 
function (x, pars = NA, fixed = FALSE, ...) 
{
    x <- as.data.frame(x)
    if (!anyNA(pars)) {
        pars <- extract_pars(pars, all_pars = names(x), fixed = fixed, 
            ...)
        x <- x[, pars, drop = FALSE]
    }
    if (!ncol(x)) {
        x <- NULL
    }
    x
}
