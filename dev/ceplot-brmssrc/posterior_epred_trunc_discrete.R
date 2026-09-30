posterior_epred_trunc_discrete <- 
function (dist, args, lb, ub) 
{
    stopifnot(is.matrix(lb), is.matrix(ub))
    message("Computing posterior_epred values for truncated ", 
        "discrete models may take a while.")
    pdf <- get(paste0("d", dist), mode = "function")
    cdf <- get(paste0("p", dist), mode = "function")
    mean_kernel <- function(x, args) {
        x * do_call(pdf, c(x, args))
    }
    if (any(is.infinite(c(lb, ub)))) {
        stop("lb and ub must be finite")
    }
    vec_lb <- lb[1, ]
    vec_ub <- ub[1, ]
    min_lb <- min(vec_lb)
    mk <- lapply((min_lb + 1):max(vec_ub), mean_kernel, args = args)
    mk <- do_call(abind, c(mk, along = 3))
    m1 <- vector("list", ncol(mk))
    for (n in seq_along(m1)) {
        J <- (vec_lb[n] - min_lb + 1):(vec_ub[n] - min_lb)
        m1[[n]] <- rowSums(mk[, n, ][, J, drop = FALSE])
    }
    rm(mk)
    m1 <- do.call(cbind, m1)
    m1/(do.call(cdf, c(list(ub), args)) - do.call(cdf, c(list(lb), 
        args)))
}
