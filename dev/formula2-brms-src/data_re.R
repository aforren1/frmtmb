function (bframe, data) 
{
    stopifnot(is.bframel(bframe))
    out <- list()
    px <- check_prefix(bframe)
    reframe <- subset2(bframe$frame$re, type = "sp", fun = "%notin%")
    if (!has_rows(reframe)) {
        return(out)
    }
    gn <- unique(reframe$gn)
    for (i in seq_along(gn)) {
        r <- subset2(reframe, gn = gn[i])
        Z <- get_model_matrix(r$form[[1]], data = data, rename = FALSE)
        idp <- paste0(r$id[1], usc(combine_prefix(px)))
        Znames <- paste0("Z_", idp, "_", r$cn)
        if (r$gtype[1] == "mm") {
            ng <- length(r$gcall[[1]]$groups)
            if (r$type[1] == "cs") {
                stop2("'cs' is not supported in multi-membership terms.")
            }
            if (r$type[1] == "mmc") {
                mmc_expr <- "^mmc\\([^:]*\\)"
                mmc_terms <- get_matches_expr(mmc_expr, colnames(Z))
                for (t in mmc_terms) {
                  pos <- which(grepl_expr(escape_all(t), colnames(Z)))
                  if (length(pos) != ng) {
                    stop2("Invalid term '", t, "': Expected ", 
                      ng, " coefficients but found ", length(pos), 
                      ".")
                  }
                  for (j in seq_along(Znames)) {
                    for (k in seq_len(ng)) {
                      out[[paste0(Znames[j], "_", k)]] <- as.array(Z[, 
                        pos[k]])
                    }
                  }
                }
            }
            else {
                for (j in seq_along(Znames)) {
                  out[paste0(Znames[j], "_", seq_len(ng))] <- list(as.array(Z[, 
                    j]))
                }
            }
        }
        else {
            if (r$type[1] == "cs") {
                ncatM1 <- nrow(r)/ncol(Z)
                Z_temp <- vector("list", ncol(Z))
                for (k in seq_along(Z_temp)) {
                  Z_temp[[k]] <- replicate(ncatM1, Z[, k], simplify = FALSE)
                }
                Z <- do_call(cbind, unlist(Z_temp, recursive = FALSE))
            }
            if (r$type[1] == "mmc") {
                stop2("'mmc' is only supported in multi-membership terms.")
            }
            for (j in seq_cols(Z)) {
                out[[Znames[j]]] <- as.array(Z[, j])
            }
        }
    }
    out
}
