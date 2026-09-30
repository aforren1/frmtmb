expand_dot_formula <- 
function (formula, data = NULL) 
{
    if (isTRUE("." %in% all.vars(formula))) {
        att <- attributes(formula)
        try_terms <- try(stats::terms(formula, data = data), 
            silent = TRUE)
        if (!is_try_error(try_terms)) {
            formula <- formula(try_terms)
        }
        attributes(formula) <- att
    }
    formula
}
