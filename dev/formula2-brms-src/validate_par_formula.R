function (formula, par = NULL, rsv_pars = NULL) 
{
    stopifnot(length(par) <= 1L)
    try_formula <- try(as_formula(formula), silent = TRUE)
    if (is_try_error(try_formula)) {
        if (length(formula) != 1L) {
            stop2("Expecting a single value when fixing parameter '", 
                par, "'.")
        }
        scalar <- SW(as.numeric(formula))
        if (!is.na(scalar)) {
            formula <- scalar
        }
        else {
            formula <- as.character(formula)
        }
        out <- named_list(par, formula)
    }
    else {
        formula <- try_formula
        if (!is.null(lhs(formula))) {
            resp_pars <- all.vars(formula[[2]])
            out <- named_list(resp_pars, list(formula))
            for (i in seq_along(out)) {
                out[[i]][[2]] <- eval2(paste("quote(", resp_pars[i], 
                  ")"))
            }
        }
        else {
            if (!isTRUE(nzchar(par))) {
                stop2("Additional formulas must be named.")
            }
            formula <- formula(paste(par, formula2str(formula)))
            out <- named_list(par, list(formula))
        }
    }
    pars <- names(out)
    if (any(grepl("\\.|_", pars))) {
        stop2("Parameter names should not contain dots or underscores.")
    }
    inv_pars <- intersect(pars, rsv_pars)
    if (length(inv_pars)) {
        stop2("The following parameter names are reserved", "for this model:\n", 
            collapse_comma(inv_pars))
    }
    out
}
