function (bterms, data, old_levels = NULL) 
{
    data <- combine_groups(data, get_group_vars(bterms))
    re <- get_re(bterms)
    out <- vector("list", nrow(re))
    used_ids <- new_ids <- NULL
    id_groups <- list()
    j <- 1
    for (i in seq_rows(re)) {
        if (!nzchar(re$type[i])) {
            coef <- colnames(get_model_matrix(re$form[[i]], data))
        }
        else if (re$type[i] == "sp") {
            coef <- frame_sp(re$form[[i]], data)$coef
        }
        else if (re$type[i] == "mmc") {
            coef <- rename(all_terms(re$form[[i]]))
        }
        else if (re$type[i] == "cs") {
            resp <- re$resp[i]
            if (nzchar(resp)) {
                stopifnot(is.mvbrmsterms(bterms))
                nthres <- max(get_thres(bterms$terms[[resp]]))
            }
            else {
                stopifnot(is.brmsterms(bterms))
                nthres <- max(get_thres(bterms))
            }
            indices <- paste0("[", seq_len(nthres), "]")
            coef <- colnames(get_model_matrix(re$form[[i]], data = data))
            coef <- as.vector(t(outer(coef, indices, paste0)))
        }
        avoid_dpars(coef, bterms)
        rdat <- data.frame(id = re$id[[i]], group = re$group[[i]], 
            gn = re$gn[[i]], gtype = re$gtype[[i]], coef = coef, 
            cn = NA, resp = re$resp[[i]], dpar = re$dpar[[i]], 
            nlpar = re$nlpar[[i]], ggn = NA, cor = re$cor[[i]], 
            type = re$type[[i]], by = re$gcall[[i]]$by, cov = re$gcall[[i]]$cov, 
            dist = re$gcall[[i]]$dist, stringsAsFactors = FALSE)
        bylevels <- NULL
        if (nzchar(rdat$by[1])) {
            bylevels <- eval2(rdat$by[1], data)
            bylevels <- rm_wsp(extract_levels(bylevels))
        }
        rdat$bylevels <- repl(bylevels, nrow(rdat))
        rdat$form <- repl(re$form[[i]], nrow(rdat))
        rdat$gcall <- repl(re$gcall[[i]], nrow(rdat))
        id <- re$id[[i]]
        if (is.na(id)) {
            rdat$id <- j
            j <- j + 1
        }
        else {
            if (id %in% used_ids) {
                k <- match(id, used_ids)
                rdat$id <- new_ids[k]
                new_id_groups <- c(re$group[[i]], re$gcall[[i]]$groups)
                if (!identical(new_id_groups, id_groups[[k]])) {
                  stop2("Can only combine group-level terms of the ", 
                    "same grouping factors.")
                }
            }
            else {
                used_ids <- c(used_ids, id)
                k <- length(used_ids)
                rdat$id <- new_ids[k] <- j
                id_groups[[k]] <- c(re$group[[i]], re$gcall[[i]]$groups)
                j <- j + 1
            }
        }
        out[[i]] <- rdat
    }
    out <- do_call(rbind, c(list(empty_reframe()), out))
    rsv_groups <- out[nzchar(out$gtype), "group"]
    other_groups <- out[!nzchar(out$gtype), "group"]
    inv_groups <- intersect(rsv_groups, other_groups)
    if (length(inv_groups)) {
        inv_groups <- paste0("'", inv_groups, "'", collapse = ", ")
        stop2("Grouping factor names ", inv_groups, " are resevered.")
    }
    dup <- duplicated(out[, c("group", "coef", vars_prefix())])
    if (any(dup)) {
        dr <- out[which(dup)[1], ]
        stop2("Duplicated group-level effects are not allowed.\n", 
            "Occured for effect '", dr$coef, "' of group '", 
            dr$group, "'.")
    }
    if (has_rows(out)) {
        for (id in unique(out$id)) {
            out$cn[out$id == id] <- seq_len(sum(out$id == id))
        }
        out$ggn <- match(out$group, unique(out$group))
        rsub <- out[!duplicated(out$group), ]
        levels <- named_list(rsub$group)
        for (i in seq_along(levels)) {
            levels[[i]] <- unique(ulapply(rsub$gcall[[i]]$groups, 
                function(g) extract_levels(get(g, data))))
            bysel <- out$group == names(levels)[i] & nzchar(out$by) & 
                !duplicated(out$by)
            bysel <- which(bysel)
            if (length(bysel) > 1L) {
                stop2("Each grouping factor can only be associated with one 'by' variable.")
            }
            if (length(bysel) == 1L) {
                rsub[i, ] <- out[bysel, ]
            }
            if (nzchar(rsub$by[i])) {
                stopifnot(rsub$type[i] %in% c("", "mmc"))
                by <- rsub$by[i]
                bylevels <- rsub$bylevels[[i]]
                byvar <- rm_wsp(eval2(by, data))
                groups <- rsub$gcall[[i]]$groups
                if (rsub$gtype[i] == "mm") {
                  byvar <- as.matrix(byvar)
                  if (!identical(dim(byvar), c(nrow(data), length(groups)))) {
                    stop2("Grouping structure 'mm' expects 'by' to be ", 
                      "a matrix with as many columns as grouping factors.")
                  }
                  df <- J <- named_list(groups)
                  for (k in seq_along(groups)) {
                    J[[k]] <- match(get(groups[k], data), levels[[i]])
                    df[[k]] <- data.frame(J = J[[k]], by = byvar[, 
                      k])
                  }
                  J <- unlist(J)
                  df <- do_call(rbind, df)
                }
                else {
                  J <- match(get(groups, data), levels[[i]])
                  df <- data.frame(J = J, by = byvar)
                }
                df <- unique(df)
                if (nrow(df) > length(unique(J))) {
                  stop2("Some levels of ", collapse_comma(groups), 
                    " correspond to multiple levels of '", by, 
                    "'.")
                }
                df <- df[order(df$J), ]
                by_per_level <- bylevels[match(df$by, bylevels)]
                attr(levels[[i]], "by") <- by_per_level
            }
        }
        if (!is.null(old_levels)) {
            set_levels(out) <- old_levels
            set_levels(out, "used") <- levels
        }
        else {
            set_levels(out) <- levels
        }
        out <- update_ranef_cov(out, bterms)
    }
    out <- out[order(out$id), , drop = FALSE]
    class(out) <- reframe_class()
    out
}
