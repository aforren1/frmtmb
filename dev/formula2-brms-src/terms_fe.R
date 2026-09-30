function (formula) 
{
    if (!is.terms(formula)) {
        formula <- terms(formula)
    }
    all_terms <- all_terms(formula)
    sp_terms <- find_terms(all_terms, "all", complete = FALSE)
    re_terms <- all_terms[grepl("\\|", all_terms)]
    int_term <- attr(formula, "intercept")
    fe_terms <- setdiff(all_terms, c(sp_terms, re_terms))
    out <- paste(c(int_term, fe_terms), collapse = "+")
    out <- str2formula(out)
    attr(out, "allvars") <- allvars_formula(out)
    attr(out, "decomp") <- get_decomp(formula)
    if (has_rsv_intercept(out, has_intercept(formula))) {
        attr(out, "int") <- FALSE
    }
    if (no_cmc(formula)) {
        attr(out, "cmc") <- FALSE
    }
    if (no_center(formula)) {
        attr(out, "center") <- FALSE
    }
    if (is_sparse(formula)) {
        attr(out, "sparse") <- TRUE
    }
    out
}
