posterior_samples.brmsfit <- 
function (x, pars = NA, fixed = FALSE, add_chain = FALSE, subset = NULL, 
    as.matrix = FALSE, as.array = FALSE, ...) 
{
    if (as.matrix && as.array) {
        stop2("Cannot use 'as.matrix' and 'as.array' at the same time.")
    }
    if (add_chain && as.array) {
        stop2("Cannot use 'add_chain' and 'as.array' at the same time.")
    }
    contains_draws(x)
    pars <- extract_pars(pars, variables(x), fixed = fixed, ...)
    iter <- x$fit@sim$iter
    warmup <- x$fit@sim$warmup
    thin <- x$fit@sim$thin
    chains <- x$fit@sim$chains
    final_iter <- ceiling((iter - warmup)/thin)
    samples_taken <- seq(warmup + 1, iter, thin)
    samples <- NULL
    if (length(pars)) {
        if (as.matrix) {
            samples <- as.matrix(x$fit, pars = pars)
        }
        else if (as.array) {
            samples <- as.array(x$fit, pars = pars)
        }
        else {
            samples <- as.data.frame(x$fit, pars = pars)
        }
        if (add_chain) {
            samples <- cbind(samples, chain = factor(rep(1:chains, 
                each = final_iter)), iter = rep(samples_taken, 
                chains))
        }
        if (!is.null(subset)) {
            if (as.array) {
                samples <- samples[subset, , , drop = FALSE]
            }
            else {
                samples <- samples[subset, , drop = FALSE]
            }
        }
    }
    samples
}
