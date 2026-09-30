prepare_predictions_fe <- 
function (bframe, draws, sdata, ...) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    if (is.null(bframe[["fe"]])) {
        return(out)
    }
    p <- usc(combine_prefix(bframe))
    X <- sdata[[paste0("X", p)]]
    fixef <- bframe$frame$fe$vars
    if (length(fixef)) {
        out$X <- X
        b_pars <- paste0("b", p, "_", fixef)
        out$b <- prepare_draws(draws, b_pars, scalar = TRUE)
    }
    out
}
