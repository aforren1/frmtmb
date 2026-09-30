function (Z, gf, max_level = NULL, weights = NULL) 
{
    if (!is.list(Z)) {
        Z <- list(Z)
    }
    if (!is.list(gf)) {
        gf <- list(gf)
    }
    if (is.null(weights)) {
        weights <- rep(1, length(gf[[1]]))
    }
    if (!is.list(weights)) {
        weights <- list(weights)
    }
    if (is.null(max_level)) {
        max_level <- max(unlist(gf))
    }
    levels <- unique(unlist(gf))
    nranef <- ncol(Z[[1]])
    Z <- mapply(expand_matrix, A = Z, x = gf, weights = weights, 
        MoreArgs = nlist(max_level))
    Z <- Reduce("+", Z)
    subset_levels(Z, levels, nranef)
}
