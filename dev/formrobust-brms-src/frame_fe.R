frame_fe <- 
function (x, data = NULL, ...) 
{
    stopifnot(is.btl(x))
    sdata <- x$sdata$fe
    if (is.null(sdata)) {
        sdata <- data_fe(x, data)
    }
    out <- list(vars = colnames(x$sdata$fe$X), center = stan_center_X(x), 
        sparse = is_sparse(x$fe), decomp = get_decomp(x$fe))
    out$vars_stan <- out$vars
    if (out$center) {
        out$vars_stan <- setdiff(out$vars_stan, "Intercept")
    }
    out
}
