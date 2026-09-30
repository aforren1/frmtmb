prepare_predictions_re_global <- 
function (bframe, draws, sdata, old_reframe, resp = NULL, sample_new_levels = "uncertainty", 
    ...) 
{
    reframe <- bframe$frame$re
    if (!has_rows(reframe)) {
        return(list())
    }
    stopifnot(is.reframe(reframe))
    resp <- resp %||% ""
    groups <- unique(reframe$group)
    old_levels <- get_levels(old_reframe)
    used_levels <- get_levels(sdata, prefix = "used")
    out <- named_list(groups, list())
    for (g in groups) {
        reframe_g <- subset2(reframe, group = g)
        old_reframe_g <- subset2(old_reframe, group = g)
        used_levels_g <- used_levels[[g]]
        old_levels_g <- old_levels[[g]]
        nlevels <- length(old_levels_g)
        nranef <- nrow(reframe_g)
        rpars <- paste0("^r_", g, "(__.+)?\\[")
        rdraws <- prepare_draws(draws, rpars, regex = TRUE)
        if (!length(rdraws)) {
            stop2("Group-level coefficients of group '", g, "' not found. ", 
                "You can control saving those coefficients via 'save_pars()'.")
        }
        cols_match <- c("coef", "resp", "dpar", "nlpar")
        used_rpars <- which(find_rows(old_reframe_g, ls = reframe_g[cols_match]))
        used_rpars <- outer(seq_len(nlevels), (used_rpars - 1) * 
            nlevels, "+")
        used_rpars <- as.vector(used_rpars)
        rdraws <- rdraws[, used_rpars, drop = FALSE]
        rdraws <- column_to_row_major_order(rdraws, nranef)
        gtype <- reframe_g$gtype[1]
        resp_g <- intersect(reframe_g$resp, resp)[1]
        id <- subset2(reframe_g, resp = resp)$id[1]
        idresp <- paste0(id, usc(resp_g))
        if (gtype == "mm") {
            ngf <- length(reframe_g$gcall[[1]]$groups)
            gf <- sdata[paste0("J_", idresp, "_", seq_len(ngf))]
            weights <- sdata[paste0("W_", idresp, "_", seq_len(ngf))]
        }
        else {
            gf <- sdata[paste0("J_", idresp)]
            weights <- list(rep(1, length(gf[[1]])))
        }
        args_new_rdraws <- nlist(reframe = reframe_g, gf, used_levels = used_levels_g, 
            old_levels = old_levels_g, rdraws = rdraws, draws, 
            sample_new_levels)
        new_rdraws <- do_call(get_new_rdraws, args_new_rdraws)
        max_level <- attr(new_rdraws, "max_level")
        gf <- attr(new_rdraws, "gf")
        rdraws <- cbind(rdraws, new_rdraws)
        levels <- unique(unlist(gf))
        rdraws <- subset_levels(rdraws, levels, nranef)
        out[[g]]$reframe <- reframe_g
        out[[g]]$rdraws <- rdraws
        out[[g]]$levels <- levels
        out[[g]]$nranef <- nranef
        out[[g]]$max_level <- max_level
        out[[g]]$gf <- gf
        out[[g]]$weights <- weights
    }
    out
}
