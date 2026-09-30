validate_newdata2 <- 
function (newdata2, object, ...) 
{
    stopifnot(is.brmsfit(object))
    bterms <- brmsterms(object$formula)
    validate_data2(newdata2, bterms = bterms, ...)
}
