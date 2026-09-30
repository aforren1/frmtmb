function (formula, check_response = TRUE, resp_rhs_all = TRUE, 
    ...) 
{
    x <- validate_formula(formula)
    mv <- isTRUE(x$mv)
    rescor <- mv && isTRUE(x$rescor)
    mecor <- isTRUE(x$mecor)
    formula <- x$formula
    family <- x$family
    y <- nlist(formula, family, mv, rescor, mecor)
    y$cov_ranef <- x$cov_ranef
    class(y) <- "brmsterms"
    y$resp <- ""
    if (check_response) {
        y$respform <- validate_resp_formula(formula, empty_ok = FALSE)
        if (mv) {
            y$resp <- terms_resp(y$respform)
        }
    }
    adforms <- terms_ad(formula, family, check_response)
    advars <- str2formula(ulapply(adforms, all_vars))
    y$adforms[names(adforms)] <- adforms
    if (!is.null(get_ad_expr(y, "thres", "gr"))) {
        attr(formula, "center") <- FALSE
        dp_classes <- dpar_class(names(x$pforms))
        mu_names <- names(x$pforms)[dp_classes == "mu"]
        for (dp in mu_names) {
            attr(x$pforms[[dp]], "center") <- FALSE
        }
    }
    if (is.mixfamily(family)) {
        mu_dpars <- paste0("mu", seq_along(family$mix))
        for (dp in mu_dpars) {
            x$pforms[[dp]] <- combine_formulas(formula, x$pforms[[dp]], 
                dp)
        }
        x$pforms <- move2start(x$pforms, mu_dpars)
        for (i in seq_along(family$mix)) {
            y$family$mix[[i]]$mix <- i
        }
    }
    else if (conv_cats_dpars(x$family)) {
        mu_dpars <- str_subset(x$family$dpars, "^mu")
        for (dp in mu_dpars) {
            x$pforms[[dp]] <- combine_formulas(formula, x$pforms[[dp]], 
                dp)
        }
        x$pforms <- move2start(x$pforms, mu_dpars)
    }
    else {
        x$pforms[["mu"]] <- combine_formulas(formula, x$pforms[["mu"]], 
            "mu")
        x$pforms <- move2start(x$pforms, "mu")
    }
    dpars <- intersect(names(x$pforms), valid_dpars(family))
    dpar_forms <- x$pforms[dpars]
    nlpars <- setdiff(names(x$pforms), dpars)
    y$dpars <- named_list(dpars)
    for (dp in dpars) {
        if (get_nl(dpar_forms[[dp]])) {
            y$dpars[[dp]] <- terms_nlf(dpar_forms[[dp]], nlpars, 
                y$resp)
        }
        else {
            y$dpars[[dp]] <- terms_lf(dpar_forms[[dp]])
        }
        y$dpars[[dp]]$family <- dpar_family(family, dp)
        y$dpars[[dp]]$dpar <- dp
        y$dpars[[dp]]$resp <- y$resp
        if (dpar_class(dp) == "mu") {
            y$dpars[[dp]]$respform <- y$respform
            y$dpars[[dp]]$adforms <- y$adforms
        }
        y$dpars[[dp]]$transform <- stan_eta_transform(y, y$dpars[[dp]]$family)
        check_cs(y$dpars[[dp]])
    }
    y$nlpars <- named_list(nlpars)
    if (length(nlpars)) {
        nlpar_forms <- x$pforms[nlpars]
        for (nlp in nlpars) {
            if (is.null(attr(nlpar_forms[[nlp]], "center"))) {
                attr(nlpar_forms[[nlp]], "center") <- FALSE
            }
            if (get_nl(nlpar_forms[[nlp]])) {
                y$nlpars[[nlp]] <- terms_nlf(nlpar_forms[[nlp]], 
                  nlpars, y$resp)
            }
            else {
                y$nlpars[[nlp]] <- terms_lf(nlpar_forms[[nlp]])
            }
            y$nlpars[[nlp]]$nlpar <- nlp
            y$nlpars[[nlp]]$resp <- y$resp
            check_cs(y$nlpars[[nlp]])
        }
        used_nlpars <- ufrom_list(c(y$dpars, y$nlpars), "used_nlpars")
        unused_nlpars <- setdiff(nlpars, used_nlpars)
        if (length(unused_nlpars)) {
            stop2("The parameter '", unused_nlpars[1], "' is not a ", 
                "valid distributional or non-linear parameter. ", 
                "Did you forget to set 'nl = TRUE'?")
        }
        used_nlpars <- from_list(y$nlpars, "used_nlpars")
        sorted_nlpars <- sort_dependencies(used_nlpars)
        y$nlpars <- y$nlpars[sorted_nlpars]
    }
    valid_dpars <- valid_dpars(y)
    inv_fixed_dpars <- setdiff(names(x$pfix), valid_dpars)
    if (length(inv_fixed_dpars)) {
        stop2("Invalid fixed parameters: ", collapse_comma(inv_fixed_dpars))
    }
    if ("sigma" %in% valid_dpars && no_sigma(y)) {
        if ("sigma" %in% c(names(x$pforms), names(x$pfix))) {
            stop2("Cannot predict or fix 'sigma' in this model.")
        }
        x$pfix$sigma <- 0
    }
    if ("nu" %in% valid_dpars && no_nu(y)) {
        if ("nu" %in% c(names(x$pforms), names(x$pfix))) {
            stop2("Cannot predict or fix 'nu' in this model.")
        }
        x$pfix$nu <- 1
    }
    disc_pars <- valid_dpars[dpar_class(valid_dpars) %in% "disc"]
    for (dp in disc_pars) {
        if (!dp %in% c(names(x$pforms), names(x$pfix))) {
            x$pfix[[dp]] <- 1
        }
    }
    for (dp in names(x$pfix)) {
        y$fdpars[[dp]] <- list(value = x$pfix[[dp]], dpar = dp)
    }
    check_fdpars(y$fdpars)
    y$unused <- attr(x$formula, "unused")
    lhsvars <- if (resp_rhs_all) 
        all_vars(y$respform)
    y$allvars <- allvars_formula(lhsvars, advars, lapply(y$dpars, 
        get_allvars), lapply(y$nlpars, get_allvars), y$time$allvars, 
        get_unused_arg_vars(y), .env = environment(formula))
    if (check_response) {
        formula_allvars <- y$respform
        formula_allvars[[3]] <- y$allvars[[2]]
        environment(formula_allvars) <- environment(y$allvars)
        y$allvars <- formula_allvars
    }
    y
}
