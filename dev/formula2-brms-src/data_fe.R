function (bframe, data) 
{
    stopifnot(is.btl(bframe))
    if (!is.null(bframe$sdata$fe)) {
        return(bframe$sdata$fe)
    }
    out <- list()
    p <- usc(combine_prefix(bframe))
    is_ord <- is_ordinal(bframe)
    cols2remove <- if (is_ord) 
        "(Intercept)"
    X <- get_model_matrix(rhs(bframe$fe), data, cols2remove = cols2remove)
    avoid_dpars(colnames(X), bframe)
    out[[paste0("K", p)]] <- ncol(X)
    if (stan_center_X(bframe)) {
        out[[paste0("Kc", p)]] <- ncol(X) - ifelse(is_ord, 0, 
            1)
    }
    out[[paste0("X", p)]] <- X
    out
}
