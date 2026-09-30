data_gr_global <- 
function (bframe, data2) 
{
    stopifnot(is.anybrmsframe(bframe))
    out <- list()
    reframe <- bframe$frame$re
    for (id in unique(reframe$id)) {
        tmp <- list()
        id_reframe <- subset2(reframe, id = id)
        nranef <- nrow(id_reframe)
        group <- id_reframe$group[1]
        levels <- attr(reframe, "levels")[[group]]
        tmp$N <- length(levels)
        tmp$M <- nranef
        tmp$NC <- as.integer(nranef * (nranef - 1)/2)
        if (nzchar(id_reframe$by[1])) {
            stopifnot(!nzchar(id_reframe$type[1]))
            bylevels <- id_reframe$bylevels[[1]]
            Jby <- match(attr(levels, "by"), bylevels)
            tmp$Nby <- length(bylevels)
            tmp$Jby <- as.array(Jby)
        }
        cov <- id_reframe$cov[1]
        if (nzchar(cov)) {
            cov_mat <- validate_recov_matrix(data2[[cov]])
            found_levels <- rownames(cov_mat)
            found <- levels %in% found_levels
            if (any(!found)) {
                stop2("Levels of the within-group covariance matrix for '", 
                  group, "' do not match names of the grouping levels.")
            }
            cov_mat <- cov_mat[levels, levels, drop = FALSE]
            tmp$Lcov <- t(chol(cov_mat))
        }
        names(tmp) <- paste0(names(tmp), "_", id)
        c(out) <- tmp
    }
    out
}
