fix_factor_contrasts <- 
function (data, olddata = NULL, ignore = NULL) 
{
    stopifnot(is(data, "data.frame"))
    stopifnot(is.null(olddata) || is.list(olddata))
    olddata <- as.data.frame(olddata)
    for (i in seq_along(data)) {
        needs_contrast <- is.factor(data[[i]]) && !names(data)[i] %in% 
            ignore
        if (needs_contrast && is.null(attr(data[[i]], "contrasts"))) {
            old_contrasts <- attr(olddata[[names(data)[i]]], 
                "contrasts")
            if (!is.null(old_contrasts)) {
                contrasts(data[[i]]) <- old_contrasts
            }
            else if (length(unique(data[[i]])) > 1) {
                contrasts(data[[i]]) <- contrasts(data[[i]])
            }
        }
    }
    data
}
