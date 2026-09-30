function (formula, data = environment(formula), cols2remove = NULL, 
    rename = TRUE, ...) 
{
    stopifnot(is_atomic_or_null(cols2remove))
    terms <- validate_terms(formula)
    if (is.null(terms)) {
        return(NULL)
    }
    if (no_int(terms)) {
        cols2remove <- union(cols2remove, "(Intercept)")
    }
    X <- stats::model.matrix(terms, data, ...)
    cols2remove <- which(colnames(X) %in% cols2remove)
    if (length(cols2remove)) {
        X <- X[, -cols2remove, drop = FALSE]
    }
    if (rename) {
        colnames(X) <- rename(colnames(X), check_dup = TRUE)
    }
    X
}
