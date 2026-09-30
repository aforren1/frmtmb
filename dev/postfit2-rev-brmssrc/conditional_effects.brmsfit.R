function (x, effects = NULL, conditions = NULL, int_conditions = NULL, 
    re_formula = NA, prob = 0.95, robust = TRUE, method = "posterior_epred", 
    spaghetti = FALSE, surface = FALSE, categorical = FALSE, 
    ordinal = FALSE, transform = NULL, resolution = 100, select_points = 0, 
    too_far = 0, probs = NULL, ...) 
{
    probs <- validate_ci_bounds(prob, probs = probs)
    method <- validate_pp_method(method)
    spaghetti <- as_one_logical(spaghetti)
    surface <- as_one_logical(surface)
    categorical <- as_one_logical(categorical)
    ordinal <- as_one_logical(ordinal)
    contains_draws(x)
    x <- restructure(x)
    new_formula <- update_re_terms(x$formula, re_formula = re_formula)
    bterms <- brmsterms(new_formula)
    if (!is.null(transform) && method != "posterior_predict") {
        stop2("'transform' is only allowed if 'method = posterior_predict'.")
    }
    if (ordinal) {
        warning2("Argument 'ordinal' is deprecated. ", "Please use 'categorical' instead.")
    }
    rsv_vars <- rsv_vars(bterms)
    use_def_effects <- is.null(effects)
    if (use_def_effects) {
        effects <- get_all_effects(bterms, rsv_vars = rsv_vars)
    }
    else {
        effects <- strsplit(as.character(effects), split = ":")
        if (any(unique(unlist(effects)) %in% rsv_vars)) {
            stop2("Variables ", collapse_comma(rsv_vars), " should not be used as effects for this model")
        }
        if (any(lengths(effects) > 2L)) {
            stop2("To display interactions of order higher than 2 ", 
                "please use the 'conditions' argument.")
        }
        all_effects <- get_all_effects(bterms, rsv_vars = rsv_vars, 
            comb_all = TRUE)
        ae_coll <- all_effects[lengths(all_effects) == 1L]
        ae_coll <- ulapply(ae_coll, paste, collapse = ":")
        matches <- match(lapply(all_effects, sort), lapply(effects, 
            sort), 0L)
        if (sum(matches) > 0 && sum(matches > 0) < length(effects)) {
            invalid <- effects[setdiff(seq_along(effects), sort(matches))]
            invalid <- ulapply(invalid, paste, collapse = ":")
            warning2("Some specified effects are invalid for this model: ", 
                collapse_comma(invalid), "\nValid effects are ", 
                "(combinations of): ", collapse_comma(ae_coll))
        }
        effects <- unique(effects[sort(matches)])
        if (!length(effects)) {
            stop2("All specified effects are invalid for this model.\n", 
                "Valid effects are (combinations of): ", collapse_comma(ae_coll))
        }
    }
    if (categorical || ordinal) {
        int_effs <- lengths(effects) == 2L
        if (any(int_effs)) {
            effects <- effects[!int_effs]
            warning2("Interactions cannot be plotted directly if 'categorical' ", 
                "is TRUE. Please use argument 'conditions' instead.")
        }
    }
    if (!length(effects)) {
        stop2("No valid effects detected.")
    }
    mf <- model.frame(x)
    conditions <- prepare_conditions(x, conditions = conditions, 
        effects = effects, re_formula = re_formula, rsv_vars = rsv_vars)
    int_conditions <- lapply(int_conditions, function(x) if (is.numeric(x)) 
        sort(x, TRUE)
    else x)
    int_vars <- get_int_vars(bterms)
    group_vars <- get_group_vars(bterms)
    out <- list()
    for (i in seq_along(effects)) {
        eff <- effects[[i]]
        cond_data <- prepare_cond_data(mf[, eff, drop = FALSE], 
            conditions = conditions, int_conditions = int_conditions, 
            int_vars = int_vars, group_vars = group_vars, surface = surface, 
            resolution = resolution, reorder = use_def_effects)
        if (surface && length(eff) == 2L && too_far > 0) {
            ex_too_far <- mgcv::exclude.too.far(g1 = cond_data[[eff[1]]], 
                g2 = cond_data[[eff[2]]], d1 = mf[, eff[1]], 
                d2 = mf[, eff[2]], dist = too_far)
            cond_data <- cond_data[!ex_too_far, ]
        }
        c(out) <- conditional_effects(bterms, fit = x, cond_data = cond_data, 
            method = method, surface = surface, spaghetti = spaghetti, 
            categorical = categorical, ordinal = ordinal, re_formula = re_formula, 
            transform = transform, conditions = conditions, int_conditions = int_conditions, 
            select_points = select_points, probs = probs, robust = robust, 
            ...)
    }
    structure(out, class = "brms_conditional_effects")
}
