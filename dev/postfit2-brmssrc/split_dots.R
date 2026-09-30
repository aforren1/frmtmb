split_dots <- 
function (x, ..., model_names = NULL, other = TRUE) 
{
    other <- as_one_logical(other)
    dots <- list(x, ...)
    names <- substitute(list(x, ...), env = parent.frame())[-1]
    names <- ulapply(names, deparse0)
    if (length(names)) {
        if (!length(names(dots))) {
            names(dots) <- names
        }
        else {
            has_no_name <- !nzchar(names(dots))
            names(dots)[has_no_name] <- names[has_no_name]
        }
    }
    is_brmsfit <- unlist(lapply(dots, is.brmsfit))
    models <- dots[is_brmsfit]
    models <- validate_models(models, model_names, names(models))
    out <- dots[!is_brmsfit]
    if (other) {
        out$models <- models
    }
    else {
        if (length(out)) {
            stop2("Only model objects can be passed to '...' for this method.")
        }
        out <- models
    }
    out
}
