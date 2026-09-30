exclude_terms.brmsfit <- 
function (x, ...) 
{
    x$formula <- exclude_terms(x$formula, ...)
    x
}
