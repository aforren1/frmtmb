nsamples.brmsfit <- 
function (object, subset = NULL, incl_warmup = FALSE, ...) 
{
    warning2("'nsamples.brmsfit' is deprecated. Please use 'ndraws' instead.")
    if (!is(object$fit, "stanfit") || !length(object$fit@sim)) {
        out <- 0
    }
    else {
        ntsamples <- object$fit@sim$n_save[1]
        if (!incl_warmup) {
            ntsamples <- ntsamples - object$fit@sim$warmup2[1]
        }
        ntsamples <- ntsamples * object$fit@sim$chains
        if (length(subset)) {
            out <- length(subset)
            if (out > ntsamples || max(subset) > ntsamples) {
                stop2("Argument 'subset' is invalid.")
            }
        }
        else {
            out <- ntsamples
        }
    }
    out
}
