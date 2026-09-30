posterior_epred_hurdle_cumulative <- 
function (prep) 
{
    adjust <- ifelse(prep$family$link == "identity", 0, 1)
    ncat_max <- max(prep$data$nthres) + adjust
    nact_min <- min(prep$data$nthres) + adjust
    init_mat <- matrix(ifelse(prep$family$link == "identity", 
        NA, 0), nrow = prep$ndraws, ncol = ncat_max - nact_min)
    args <- list(link = prep$family$link)
    out <- vector("list", prep$nobs)
    for (i in seq_along(out)) {
        args_i <- args
        args_i$eta <- slice_col(get_dpar(prep, "mu", i))
        args_i$disc <- slice_col(get_dpar(prep, "disc", i))
        args_i$thres <- subset_thres(prep, i)
        ncat_i <- NCOL(args_i$thres) + adjust
        args_i$x <- seq_len(ncat_i)
        out[[i]] <- do_call(dcumulative, args_i)
        if (ncat_i < ncat_max) {
            sel <- seq_len(ncat_max - ncat_i)
            out[[i]] <- cbind(out[[i]], init_mat[, sel])
        }
        hu <- get_dpar(prep, "hu", i)
        out[[i]] <- cbind(hu, out[[i]] * (1 - hu))
    }
    out <- abind(out, along = 3)
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- c(paste0(0), seq_len(ncat_max))
    out
}
