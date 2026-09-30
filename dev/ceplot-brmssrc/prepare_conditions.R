prepare_conditions <- 
function (fit, conditions = NULL, effects = NULL, re_formula = NA, 
    rsv_vars = NULL) 
{
    mf <- model.frame(fit)
    new_formula <- update_re_terms(fit$formula, re_formula = re_formula)
    bterms <- brmsterms(new_formula)
    if (any(grepl_expr("^(as\\.)?factor(.+)$", bterms$allvars))) {
        warning2("Using 'factor' or 'as.factor' in the model formula ", 
            "might lead to problems in 'conditional_effects'.", 
            "Please convert your variables to factors beforehand.")
    }
    req_vars <- all_vars(rhs(bterms$allvars))
    req_vars <- setdiff(req_vars, rsv_vars)
    req_vars <- setdiff(req_vars, names(fit$data2))
    if (is.null(conditions)) {
        conditions <- as.data.frame(as.list(rep(NA, length(req_vars))))
        names(conditions) <- req_vars
    }
    else {
        conditions <- as.data.frame(conditions)
        if (!nrow(conditions)) {
            stop2("Argument 'conditions' must have a least one row.")
        }
        conditions <- unique(conditions)
        if (any(duplicated(get_cond__(conditions)))) {
            stop2("Condition labels should be unique.")
        }
        req_vars <- setdiff(req_vars, names(conditions))
    }
    trial_vars <- all_vars(bterms$adforms$trials)
    trial_vars <- trial_vars[!vars_specified(trial_vars, conditions)]
    if (length(trial_vars)) {
        message("Setting all 'trials' variables to 1 by ", "default if not specified otherwise.")
        req_vars <- setdiff(req_vars, trial_vars)
        for (v in trial_vars) {
            conditions[[v]] <- 1
        }
    }
    subset_vars <- get_ad_vars(bterms, "subset")
    int_vars <- get_int_vars(bterms)
    group_vars <- get_group_vars(bterms)
    req_vars <- setdiff(req_vars, group_vars)
    for (v in req_vars) {
        if (is_like_factor(mf[[v]])) {
            if (v %in% subset_vars) {
                conditions[[v]] <- TRUE
            }
            else {
                levels <- levels(as.factor(mf[[v]]))
                ordered <- is.ordered(mf[[v]])
                conditions[[v]] <- factor(levels[1], levels, 
                  ordered = ordered)
            }
        }
        else {
            if (v %in% subset_vars) {
                conditions[[v]] <- 1
            }
            else if (v %in% int_vars) {
                conditions[[v]] <- round(median(mf[[v]], na.rm = TRUE))
            }
            else {
                conditions[[v]] <- mean(mf[[v]], na.rm = TRUE)
            }
        }
    }
    all_vars <- c(all_vars(bterms$allvars), "cond__")
    unused_vars <- setdiff(names(conditions), all_vars)
    if (length(unused_vars)) {
        warning2("The following variables in 'conditions' are not ", 
            "part of the model:\n", collapse_comma(unused_vars))
    }
    cond__ <- conditions$cond__
    conditions <- validate_newdata(conditions, fit, re_formula = re_formula, 
        allow_new_levels = TRUE, check_response = FALSE, incl_autocor = FALSE)
    conditions$cond__ <- cond__
    conditions
}
