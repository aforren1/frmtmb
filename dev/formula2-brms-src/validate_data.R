function (data, bterms, data2 = list(), knots = NULL, na_action = na_omit, 
    drop_unused_levels = TRUE, attr_terms = NULL, data_name = "") 
{
    if (missing(data)) {
        stop2("Data must be specified using the 'data' argument.")
    }
    if (is.null(knots)) {
        knots <- get_knots(data)
    }
    data <- try(as.data.frame(data), silent = TRUE)
    if (is_try_error(data)) {
        stop2("Argument 'data' must be coercible to a data.frame.")
    }
    if (!isTRUE(nrow(data) > 0L)) {
        stop2("Argument 'data' does not contain observations.")
    }
    data <- data_rsv_intercept(data, bterms = bterms)
    all_vars_formula <- bterms$allvars
    missing_vars <- setdiff(all_vars(all_vars_formula), names(data))
    if (length(missing_vars)) {
        missing_vars2 <- setdiff(missing_vars, names(data2))
        if (length(missing_vars2)) {
            stop2("The following variables can neither be found in ", 
                "'data' nor in 'data2':\n", collapse_comma(missing_vars2))
        }
        missing_vars_formula <- paste0(". ~ . ", collapse(" - ", 
            missing_vars))
        all_vars_formula <- update(all_vars_formula, missing_vars_formula)
    }
    all_vars_terms <- terms(all_vars_formula)
    terms_env <- environment(all_vars_terms)
    environment(all_vars_terms) <- as.environment(as.list(data2))
    parent.env(environment(all_vars_terms)) <- terms_env
    attributes(all_vars_terms)[names(attr_terms)] <- attr_terms
    attr(data, "terms") <- NULL
    na_action_bterms <- function(object, ...) {
        na_action(object, bterms = bterms, ...)
    }
    data <- model.frame(all_vars_terms, data, na.action = na_action_bterms, 
        drop.unused.levels = drop_unused_levels)
    if (any(grepl("__|_$", colnames(data)))) {
        stop2("Variable names may not contain double underscores ", 
            "or underscores at the end.")
    }
    if (!isTRUE(nrow(data) > 0L)) {
        stop2("All observations in the data were removed ", "presumably because of NA values.")
    }
    if (any(ulapply(data, is.infinite))) {
        warning2("Found infinite values in the data, ", "which may cause issues for Stan.")
    }
    groups <- get_group_vars(bterms)
    data <- combine_groups(data, groups)
    data <- fix_factor_contrasts(data, ignore = groups)
    data <- order_data(data, bterms = bterms)
    attr(data, "knots") <- knots
    attr(data, "drop_unused_levels") <- drop_unused_levels
    attr(data, "data_name") <- data_name
    data
}
