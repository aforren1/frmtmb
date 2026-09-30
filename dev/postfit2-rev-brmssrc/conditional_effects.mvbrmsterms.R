function (x, resp = NULL, ...) 
{
    resp <- validate_resp(resp, x$responses)
    x$terms <- x$terms[resp]
    out <- lapply(x$terms, conditional_effects, ...)
    unlist(out, recursive = FALSE)
}
