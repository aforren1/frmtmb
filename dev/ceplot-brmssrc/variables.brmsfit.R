variables.brmsfit <- 
function (x, ...) 
{
    out <- dimnames(x$fit)
    if (is.list(out)) {
        out <- out$parameters
    }
    out
}
