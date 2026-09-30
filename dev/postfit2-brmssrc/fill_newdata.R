fill_newdata <- 
function (newdata, vars, olddata = NULL, n = 1) 
{
    stopifnot(is.data.frame(newdata), is.character(vars))
    vars <- setdiff(vars, names(newdata))
    if (is.null(olddata)) {
        if (length(vars)) {
            newdata[, vars] <- NA
        }
        return(newdata)
    }
    stopifnot(is.data.frame(olddata), length(n) == 1)
    for (v in vars) {
        cval <- olddata[n, v] %||% NA
        if (length(dim(cval)) == 2) {
            cval <- matrix(cval, nrow(newdata), ncol(cval), byrow = TRUE)
        }
        newdata[[v]] <- cval
    }
    newdata
}
