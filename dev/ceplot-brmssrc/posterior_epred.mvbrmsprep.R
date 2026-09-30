posterior_epred.mvbrmsprep <- 
function (object, ...) 
{
    out <- lapply(object$resps, posterior_epred, ...)
    along <- ifelse(length(out) > 1, 3, 2)
    do_call(abind, c(out, along = along))
}
