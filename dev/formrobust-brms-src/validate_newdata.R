validate_newdata <- 
function (newdata, object, re_formula = NULL, allow_new_levels = FALSE, 
    newdata2 = NULL, resp = NULL, check_response = TRUE, incl_autocor = TRUE, 
    group_vars = NULL, req_vars = NULL, ...) 
{
    newdata <- try(as.data.frame(newdata), silent = TRUE)
    if (is_try_error(newdata)) {
        stop2("Argument 'newdata' must be coercible to a data.frame.")
    }
    object <- restructure(object)
    object <- exclude_terms(object, incl_autocor = incl_autocor)
    resp <- validate_resp(resp, object)
    new_formula <- update_re_terms(formula(object), re_formula)
    bterms <- brmsterms(new_formula, resp_rhs_all = FALSE)
    all_vars <- all.vars(bterms$allvars)
    if (is.null(req_vars)) {
        req_vars <- all_vars
    }
    else {
        req_vars <- as.character(req_vars)
        req_vars <- intersect(req_vars, all_vars)
    }
    if (is.mvbrmsterms(bterms) && !is.null(resp)) {
        resp <- validate_resp(resp, bterms$responses)
        form_req_vars <- from_list(bterms$terms[resp], "allvars")
        form_req_vars <- allvars_formula(form_req_vars)
        req_vars <- intersect(req_vars, all.vars(form_req_vars))
    }
    not_req_vars <- setdiff(all_vars, req_vars)
    not_req_vars <- setdiff(not_req_vars, names(newdata))
    newdata <- fill_newdata(newdata, not_req_vars, object$data)
    only_resp <- all.vars(bterms$respform)
    only_resp <- setdiff(only_resp, all.vars(rhs(bterms$allvars)))
    dec_vars <- get_ad_vars(bterms, "dec")
    missing_resp <- setdiff(c(only_resp, dec_vars), names(newdata))
    if (length(missing_resp)) {
        if (check_response) {
            stop2("Response variables must be specified in 'newdata'.\n", 
                "Missing variables: ", collapse_comma(missing_resp))
        }
        else {
            newdata <- fill_newdata(newdata, missing_resp)
        }
    }
    cens_vars <- get_ad_vars(bterms, "cens")
    for (v in setdiff(cens_vars, names(newdata))) {
        newdata[[v]] <- 0
    }
    weights_vars <- get_ad_vars(bterms, "weights")
    for (v in setdiff(weights_vars, names(newdata))) {
        newdata[[v]] <- 1
    }
    mf <- model.frame(object)
    for (i in seq_along(mf)) {
        if (is_like_factor(mf[[i]])) {
            mf[[i]] <- as.factor(mf[[i]])
        }
    }
    pw_vars <- ufrom_list(get_re(bterms)$gcall, "pw")
    for (v in setdiff(pw_vars, names(newdata))) {
        newdata[[v]] <- 1
    }
    newdata <- data_rsv_intercept(newdata, bterms)
    new_group_vars <- get_group_vars(bterms)
    if (allow_new_levels && length(new_group_vars)) {
        mis_group_vars <- new_group_vars[!grepl(":", new_group_vars)]
        mis_group_vars <- setdiff(mis_group_vars, names(newdata))
        newdata <- fill_newdata(newdata, mis_group_vars)
    }
    newdata <- combine_groups(newdata, new_group_vars)
    if (is.null(group_vars)) {
        group_vars <- get_group_vars(object)
    }
    do_check <- union(get_pred_vars(bterms), get_int_vars(bterms))
    unused_arg_vars <- get_unused_arg_vars(bterms)
    dont_check <- unique(c(group_vars, cens_vars, unused_arg_vars))
    dont_check <- setdiff(dont_check, do_check)
    dont_check <- names(mf) %in% dont_check
    is_factor <- ulapply(mf, is.factor)
    factors <- mf[is_factor & !dont_check]
    if (length(factors)) {
        factor_names <- names(factors)
        for (i in seq_along(factors)) {
            new_factor <- newdata[[factor_names[i]]]
            if (!is.null(new_factor)) {
                if (!is.factor(new_factor)) {
                  new_factor <- factor(new_factor)
                }
                old_levels <- levels(factors[[i]])
                if (length(old_levels) <= 1) {
                  next
                }
                new_levels <- levels(new_factor)
                old_contrasts <- contrasts(factors[[i]])
                old_ordered <- is.ordered(factors[[i]])
                to_zero <- is.na(new_factor) | new_factor %in% 
                  "zero__"
                is_resp <- factor_names[i] %in% all.vars(bterms$respform)
                if (!is_resp && any(to_zero)) {
                  levels(new_factor) <- c(new_levels, "zero__")
                  new_factor[to_zero] <- "zero__"
                  old_levels <- c(old_levels, "zero__")
                  old_contrasts <- rbind(old_contrasts, zero__ = 0)
                }
                if (any(!new_levels %in% old_levels)) {
                  stop2("New factor levels are not allowed.", 
                    "\nLevels allowed: ", collapse_comma(old_levels), 
                    "\nLevels found: ", collapse_comma(new_levels))
                }
                newdata[[factor_names[i]]] <- factor(new_factor, 
                  old_levels, ordered = old_ordered)
                attr(newdata[[factor_names[i]]], "contrasts") <- old_contrasts
            }
        }
    }
    num_names <- names(mf)[!is_factor]
    num_names <- setdiff(num_names, group_vars)
    for (nm in intersect(num_names, names(newdata))) {
        if (!anyNA(newdata[[nm]]) && !is.numeric(newdata[[nm]])) {
            stop2("Variable '", nm, "' was originally ", "numeric but is not in 'newdata'.")
        }
    }
    mo_vars <- get_sp_vars(bterms, "mo")
    if (length(mo_vars)) {
        num_mo_vars <- names(mf)[!is_factor & names(mf) %in% 
            mo_vars]
        for (v in num_mo_vars) {
            new_values <- get(v, newdata)
            min_value <- min(mf[[v]])
            invalid <- new_values < min_value | new_values > 
                max(mf[[v]])
            invalid <- invalid | !is_wholenumber(new_values)
            if (sum(invalid)) {
                stop2("Invalid values in variable '", v, "': ", 
                  collapse_comma(new_values[invalid]))
            }
            attr(newdata[[v]], "min") <- min_value
        }
    }
    used_vars <- c(names(newdata), all.vars(bterms$allvars))
    used_vars <- union(used_vars, rsv_vars(bterms))
    all_vars <- all.vars(str2formula(names(mf)))
    unused_vars <- setdiff(all_vars, used_vars)
    newdata <- fill_newdata(newdata, unused_vars)
    old_levels <- get_levels(bterms, data = mf)
    if (!allow_new_levels) {
        new_levels <- get_levels(bterms, data = newdata)
        for (g in names(old_levels)) {
            unknown_levels <- setdiff(new_levels[[g]], old_levels[[g]])
            if (anyNA(newdata[[g]])) {
                c(unknown_levels) <- NA
            }
            if (length(unknown_levels)) {
                unknown_levels <- collapse_comma(unknown_levels)
                stop2("Levels ", unknown_levels, " of grouping factor '", 
                  g, "' ", "cannot be found in the fitted model. ", 
                  "Consider setting argument 'allow_new_levels' to TRUE.")
            }
        }
    }
    old_terms <- attr(object$data, "terms")
    attr_terms <- c("variables", "predvars")
    attr_terms <- attributes(old_terms)[attr_terms]
    newdata <- validate_data(newdata, bterms = bterms, na_action = na.pass, 
        drop_unused_levels = FALSE, attr_terms = attr_terms, 
        data2 = current_data2(object, newdata2), knots = get_knots(object$data))
    newdata
}
