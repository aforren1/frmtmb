get_all_effects.mvbrmsterms <- 
function (x, ...) 
{
    out <- lapply(x$terms, get_all_effects, ...)
    unique(unlist(out, recursive = FALSE))
}
