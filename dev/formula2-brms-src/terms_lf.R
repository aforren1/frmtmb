function (formula) 
{
    formula <- rhs(as.formula(formula))
    y <- nlist(formula)
    formula <- terms(formula)
    check_accidental_helper_functions(formula)
    types <- setdiff(all_term_types(), excluded_term_types(formula))
    for (t in types) {
        tmp <- do_call(paste0("terms_", t), list(formula))
        if (is.data.frame(tmp) || is.formula(tmp)) {
            y[[t]] <- tmp
        }
    }
    y$allvars <- allvars_formula(get_allvars(y$fe), get_allvars(y$re), 
        get_allvars(y$cs), get_allvars(y$sp), get_allvars(y$sm), 
        get_allvars(y$gp), get_allvars(y$ac), get_allvars(y$offset))
    structure(y, class = "btl")
}
