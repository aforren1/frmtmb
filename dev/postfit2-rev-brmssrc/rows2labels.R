function (x, digits = 2, sep = " & ", incl_vars = TRUE, ...) 
{
    x <- as.data.frame(x)
    incl_vars <- as_one_logical(incl_vars)
    out <- x
    for (i in seq_along(out)) {
        if (!is_like_factor(out[[i]])) {
            out[[i]] <- round(out[[i]], digits)
        }
        if (incl_vars) {
            out[[i]] <- paste0(names(out)[i], " = ", out[[i]])
        }
    }
    paste_sep <- function(..., sep__ = sep) {
        paste(..., sep = sep__)
    }
    Reduce(paste_sep, out)
}
