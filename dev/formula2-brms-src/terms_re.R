function (formula) 
{
    re_terms <- get_re_terms(formula, brackets = FALSE)
    if (!length(re_terms)) {
        return(NULL)
    }
    re_terms <- split_re_terms(re_terms)
    re_parts <- re_parts(re_terms)
    out <- allvars <- vector("list", length(re_terms))
    type <- attr(re_terms, "type")
    for (i in seq_along(re_terms)) {
        gcall <- eval2(re_parts$rhs[i])
        form <- str2formula(re_parts$lhs[i])
        group <- paste0(gcall$type, collapse(gcall$groups))
        out[[i]] <- data.frame(group = group, gtype = gcall$type, 
            gn = i, id = gcall$id, type = type[i], cor = gcall$cor, 
            stringsAsFactors = FALSE)
        out[[i]]$gcall <- list(gcall)
        out[[i]]$form <- list(form)
        ftype <- str_if(type[i] %in% "cs", "", type[i])
        re_allvars <- get_allvars(form, type = ftype)
        allvars[[i]] <- allvars_formula(re_allvars, gcall$allvars)
    }
    out <- do_call(rbind, out)
    out <- out[order(out$group), ]
    attr(out, "allvars") <- allvars_formula(allvars)
    if (no_cmc(formula)) {
        for (i in seq_rows(out)) {
            attr(out$form[[i]], "cmc") <- FALSE
        }
    }
    out
}
