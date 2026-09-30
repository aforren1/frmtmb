deparse0 <- 
function (x, max_char = NULL, ...) 
{
    out <- collapse(deparse(x, ...))
    if (isTRUE(max_char > 0)) {
        out <- substr(out, 1, max_char)
    }
    out
}
