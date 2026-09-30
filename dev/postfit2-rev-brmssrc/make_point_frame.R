function (bterms, mf, effects, conditions, select_points = 0, 
    transform = NULL, ...) 
{
    stopifnot(is.brmsterms(bterms), is.data.frame(mf))
    effects <- intersect(effects, names(mf))
    points <- mf[, effects, drop = FALSE]
    points$resp__ <- model.response(model.frame(bterms$respform, 
        mf, na.action = na.pass))
    req_vars <- names(mf)
    groups <- get_re_group_vars(bterms)
    if (length(groups)) {
        c(req_vars) <- unlist(strsplit(groups, ":"))
    }
    req_vars <- unique(setdiff(req_vars, effects))
    req_vars <- intersect(req_vars, names(conditions))
    if (length(req_vars)) {
        cond__ <- get_cond__(conditions)
        mf <- mf[, req_vars, drop = FALSE]
        conditions <- conditions[, req_vars, drop = FALSE]
        points$cond__ <- NA
        points <- replicate(nrow(conditions), points, simplify = FALSE)
        for (i in seq_along(points)) {
            cond <- conditions[i, , drop = FALSE]
            not_na <- function(x) !any(is.na(x) | x %in% "zero__")
            not_na <- ulapply(cond, not_na)
            cond <- cond[, not_na, drop = FALSE]
            mf_tmp <- mf[, not_na, drop = FALSE]
            if (ncol(mf_tmp)) {
                is_num <- sapply(mf_tmp, is.numeric)
                is_num <- is_num & !names(mf_tmp) %in% groups
                if (sum(is_num)) {
                  stopifnot(select_points >= 0)
                  if (select_points > 0) {
                    for (v in names(mf_tmp)[is_num]) {
                      min <- min(mf_tmp[, v], na.rm = TRUE)
                      max <- max(mf_tmp[, v], na.rm = TRUE)
                      unit <- scale_unit(mf_tmp[, v], min, max)
                      unit_cond <- scale_unit(cond[, v], min, 
                        max)
                      unit_diff <- abs(unit - unit_cond)
                      close_enough <- unit_diff <= select_points
                      mf_tmp[[v]][close_enough] <- cond[, v]
                      mf_tmp[[v]][!close_enough] <- NA
                    }
                  }
                  else {
                    cond <- cond[, !is_num, drop = FALSE]
                    mf_tmp <- mf_tmp[, !is_num, drop = FALSE]
                  }
                }
            }
            if (ncol(mf_tmp)) {
                K <- do_call("paste", c(mf_tmp, sep = "\r")) %in% 
                  do_call("paste", c(cond, sep = "\r"))
            }
            else {
                K <- seq_rows(mf)
            }
            points[[i]]$cond__[K] <- cond__[i]
        }
        points <- do_call(rbind, points)
        points <- points[!is.na(points$cond__), , drop = FALSE]
        points$cond__ <- factor(points$cond__, cond__)
    }
    points <- add_effects__(points, effects)
    if (!is.numeric(points$resp__)) {
        points$resp__ <- as.numeric(as.factor(points$resp__))
        if (is_binary(bterms$family)) {
            points$resp__ <- points$resp__ - 1
        }
    }
    if (!is.null(transform)) {
        points$resp__ <- do_call(transform, list(points$resp__))
    }
    points
}
