exclude_terms.mvbrmsformula <- 
function (x, ...) 
{
    for (i in seq_along(x$forms)) {
        x$forms[[i]] <- exclude_terms(x$forms[[i]], ...)
    }
    x
}
