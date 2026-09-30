posterior_epred_ordinal <- 
function (prep) 
{
    dens <- get(paste0("d", prep$family$family), mode = "function")
    adjust <- ifelse(prep$family$link == "identity", 0, 1)
    ncat_max <- max(prep$data$nthres) + adjust
    nact_min <- min(prep$data$nthres) + adjust
    init_mat <- matrix(ifelse(prep$family$link == "identity", 
        NA, 0), nrow = prep$ndraws, ncol = ncat_max - nact_min)
    args <- list(link = prep$family$link)
    out <- vector("list", prep$nobs)
    for (i in seq_along(out)) {
        args_i <- args
        args_i$eta <- slice_col(prep$dpars$mu, i)
        args_i$disc <- slice_col(prep$dpars$disc, i)
        args_i$thres <- subset_thres(prep, i)
        ncat_i <- NCOL(args_i$thres) + adjust
        args_i$x <- seq_len(ncat_i)
        out[[i]] <- do_call(dens, args_i)
        if (ncat_i < ncat_max) {
            sel <- seq_len(ncat_max - ncat_i)
            out[[i]] <- cbind(out[[i]], init_mat[, sel])
        }
    }
    out <- abind(out, along = 3)
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- seq_len(ncat_max)
    out
}
