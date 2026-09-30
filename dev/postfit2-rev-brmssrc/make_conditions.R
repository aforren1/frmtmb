function (x, vars, ...) 
{
    vars <- rev(as.character(vars))
    if (!is.data.frame(x) && "data" %in% names(x)) {
        x <- x$data
    }
    x <- as.data.frame(x)
    out <- named_list(vars)
    for (v in vars) {
        tmp <- get(v, x)
        if (is_like_factor(tmp)) {
            tmp <- levels(as.factor(tmp))
        }
        else {
            tmp <- mean(tmp, na.rm = TRUE) + (-1:1) * sd(tmp, 
                na.rm = TRUE)
        }
        out[[v]] <- tmp
    }
    out <- rev(expand.grid(out))
    out$cond__ <- rows2labels(out, ...)
    out
}
