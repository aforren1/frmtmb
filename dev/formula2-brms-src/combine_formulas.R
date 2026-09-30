function (formula1, formula2, lhs = "", update = FALSE) 
{
    stopifnot(is.formula(formula1))
    stopifnot(is.null(formula2) || is.formula(formula2))
    lhs <- as_one_character(lhs)
    update <- as_one_logical(update)
    if (is.null(formula2)) {
        rhs <- str_rhs(formula1)
        att <- attributes(formula1)
    }
    else if (update && has_terms(formula1)) {
        if (get_nl(formula1) || get_nl(formula2)) {
            stop2("Cannot combine non-linear formulas.")
        }
        old_formula <- eval2(paste0("~ ", str_rhs(formula1)))
        new_formula <- eval2(paste0("~ . + ", str_rhs(formula2)))
        rhs <- str_rhs(update(old_formula, new_formula))
        att <- attributes(formula1)
        att[names(attributes(formula2))] <- attributes(formula2)
    }
    else {
        rhs <- str_rhs(formula2)
        att <- attributes(formula2)
    }
    out <- eval2(paste0(lhs, " ~ ", rhs))
    attributes(out)[names(att)] <- att
    out
}
