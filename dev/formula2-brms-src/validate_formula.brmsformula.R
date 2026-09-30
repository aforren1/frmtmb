function (formula, family = gaussian(), autocor = NULL, data = NULL, 
    threshold = NULL, sparse = NULL, cov_ranef = NULL, ...) 
{
    out <- bf(formula)
    if (is.null(out$family) && !is.null(family)) {
        out$family <- validate_family(family)
    }
    out$formula <- expand_dot_formula(out$formula, data)
    for (i in seq_along(out$pforms)) {
        out$pforms[[i]] <- expand_dot_formula(out$pforms[[i]], 
            data)
    }
    out$mecor <- default_mecor(out$mecor)
    if (has_cat(out) && !is.null(data)) {
        out$family$cats <- extract_cat_names(out, data)
    }
    if (is_ordinal(out$family)) {
        try_terms <- try(stats::terms(out$formula), silent = TRUE)
        intercept <- attr(try_terms, "intercept", TRUE)
        if (!is_try_error(try_terms) && isTRUE(intercept == 0)) {
            stop2("Cannot remove the intercept in an ordinal model.")
        }
        if (!is.null(data)) {
            out$family$thres <- extract_thres_names(out, data)
            out$family$cats <- extract_cat_names(out, data)
        }
    }
    conv_cats_dpars <- conv_cats_dpars(out$family)
    if (conv_cats_dpars && !is.null(data)) {
        if (length(out$family$cats) < 2L) {
            stop2("At least 2 response categories are required.")
        }
        if (is.null(out$family$refcat)) {
            out$family$refcat <- out$family$cats[1]
        }
        if (isNA(out$family$refcat)) {
            predcats <- out$family$cats
        }
        else {
            if (!out$family$refcat %in% out$family$cats) {
                stop2("The reference response category must be one of ", 
                  collapse_comma(out$family$cats), ".")
            }
            predcats <- setdiff(out$family$cats, out$family$refcat)
        }
        multi_dpars <- valid_dpars(out$family, type = "multi")
        for (dp in rev(multi_dpars)) {
            dp_dpars <- make_stan_names(paste0(dp, predcats))
            if (any(duplicated(dp_dpars))) {
                stop2("Invalid response category names. Please avoid ", 
                  "using any special characters in the names.")
            }
            old_dp_dpars <- str_subset(out$family$dpars, paste0("^", 
                dp))
            out$family$dpars <- setdiff(out$family$dpars, old_dp_dpars)
            out$family$dpars <- union(dp_dpars, out$family$dpars)
        }
    }
    if (is_cox(out$family) && !is.null(data)) {
        out$family$bhaz <- extract_bhaz(out, data)
    }
    if (is.mixfamily(out$family)) {
        for (i in seq_along(out$family$mix)) {
            for (term in c("cats", "thres", "bhaz")) {
                out$family$mix[[i]][[term]] <- out$family[[term]]
            }
        }
    }
    require_threshold <- is_ordinal(out$family) && is.null(out$family$threshold)
    if (require_threshold && !is.null(threshold)) {
        out$family <- validate_family(out$family, threshold = threshold)
    }
    if (!is.null(sparse)) {
        warning2("Argument 'sparse' should be specified within the ", 
            "'formula' argument. See ?brmsformula for help.")
        sparse <- as_one_logical(sparse)
        if (is.null(attr(out$formula, "sparse"))) {
            attr(out$formula, "sparse") <- sparse
        }
        for (i in seq_along(out$pforms)) {
            if (is.null(attr(out$pforms[[i]], "sparse"))) {
                attr(out$pforms[[i]], "sparse") <- sparse
            }
        }
    }
    if (is.null(attr(out$formula, "autocor")) && !is.null(autocor)) {
        warning2("Argument 'autocor' should be specified within the ", 
            "'formula' argument. See ?brmsformula for help.")
        attr(out$formula, "autocor") <- validate_autocor(autocor)
    }
    if (!is.null(cov_ranef)) {
        out$cov_ranef <- validate_cov_ranef(cov_ranef)
    }
    bf(out)
}
