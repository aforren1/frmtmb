function (bframe, sdata, prep_re = list(), sample_new_levels = "uncertainty", 
    ...) 
{
    out <- list()
    if (!length(prep_re)) {
        return(out)
    }
    px <- check_prefix(bframe)
    p <- usc(combine_prefix(px))
    reframe_px <- from_list(prep_re, "reframe")
    reframe_px <- do_call(rbind, reframe_px)
    reframe_px <- subset2(reframe_px, ls = px)
    if (!has_rows(reframe_px)) {
        return(out)
    }
    groups <- unique(reframe_px$group)
    out[c("Z", "Zsp", "Zcs")] <- list(named_list(groups))
    for (g in groups) {
        reframe_g <- prep_re[[g]]$reframe
        reframe_g_px <- subset2(reframe_g, ls = px)
        rdraws <- prep_re[[g]]$rdraws
        nranef <- prep_re[[g]]$nranef
        levels <- prep_re[[g]]$levels
        max_level <- prep_re[[g]]$max_level
        gf <- prep_re[[g]]$gf
        weights <- prep_re[[g]]$weights
        reframe_g_px_sp <- subset2(reframe_g_px, type = "sp")
        if (nrow(reframe_g_px_sp)) {
            Z <- matrix(1, length(gf[[1]]))
            out[["Zsp"]][[g]] <- prepare_Z(Z, gf, max_level, 
                weights)
            for (co in reframe_g_px_sp$coef) {
                select <- find_rows(reframe_g, ls = px) & reframe_g$coef == 
                  co & reframe_g$type == "sp"
                select <- which(select)
                select <- select + nranef * (seq_along(levels) - 
                  1)
                out[["rsp"]][[co]][[g]] <- rdraws[, select, drop = FALSE]
            }
        }
        reframe_g_px_cs <- subset2(reframe_g_px, type = "cs")
        if (nrow(reframe_g_px_cs)) {
            reframe_g_px_cs_1 <- reframe_g_px_cs[grepl("\\[1\\]$", 
                reframe_g_px_cs$coef), ]
            Znames <- paste0("Z_", reframe_g_px_cs_1$id, p, "_", 
                reframe_g_px_cs_1$cn)
            Z <- do_call(cbind, sdata[Znames])
            out[["Zcs"]][[g]] <- prepare_Z(Z, gf, max_level, 
                weights)
            for (i in seq_len(sdata$nthres)) {
                index <- paste0("\\[", i, "\\]$")
                select <- find_rows(reframe_g, ls = px) & grepl(index, 
                  reframe_g$coef) & reframe_g$type == "cs"
                select <- which(select)
                select <- as.vector(outer(select, nranef * (seq_along(levels) - 
                  1), "+"))
                out[["rcs"]][[g]][[i]] <- rdraws[, select, drop = FALSE]
            }
        }
        reframe_g_px_basic <- subset2(reframe_g_px, type = c("", 
            "mmc"))
        if (nrow(reframe_g_px_basic)) {
            Znames <- paste0("Z_", reframe_g_px_basic$id, p, 
                "_", reframe_g_px_basic$cn)
            if (reframe_g_px_basic$gtype[1] == "mm") {
                ng <- length(reframe_g_px_basic$gcall[[1]]$groups)
                Z <- vector("list", ng)
                for (k in seq_len(ng)) {
                  Z[[k]] <- do_call(cbind, sdata[paste0(Znames, 
                    "_", k)])
                }
            }
            else {
                Z <- do_call(cbind, sdata[Znames])
            }
            out[["Z"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
            select <- find_rows(reframe_g, ls = px) & reframe_g$type %in% 
                c("", "mmc")
            select <- which(select)
            select <- as.vector(outer(select, nranef * (seq_along(levels) - 
                1), "+"))
            out[["r"]][[g]] <- rdraws[, select, drop = FALSE]
        }
    }
    out
}
