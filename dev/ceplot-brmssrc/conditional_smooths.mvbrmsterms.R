conditional_smooths.mvbrmsterms <- 
function (x, ...) 
{
    out <- list()
    for (r in names(x$terms)) {
        c(out) <- conditional_smooths(x$terms[[r]], ...)
    }
    out
}
