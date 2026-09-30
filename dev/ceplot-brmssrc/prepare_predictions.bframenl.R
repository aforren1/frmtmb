prepare_predictions.bframenl <- 
function (x, draws, sdata, ...) 
{
    out <- list(family = x$family, nlform = x$formula[[2]], env = env_stan_functions(parent = environment(x$formula)), 
        ndraws = nrow(draws), nobs = sdata[[paste0("N", usc(x$resp))]], 
        used_nlpars = x$used_nlpars, loop = x$loop)
    class(out) <- "bprepnl"
    p <- usc(combine_prefix(x))
    covars <- all.vars(x$covars)
    for (i in seq_along(covars)) {
        cvalues <- sdata[[paste0("C", p, "_", i)]]
        cdim <- c(out$ndraws, out$nobs)
        if (is.matrix(cvalues)) {
            c(cdim) <- dim(cvalues)[2]
        }
        out$C[[covars[i]]] <- data2draws(cvalues, dim = cdim)
    }
    out
}
