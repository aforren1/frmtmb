function (formula, family = NULL, autocor = NULL, cov_ranef = NULL, 
    ...) 
{
    nresp <- length(formula$forms)
    if (!is(family, "list")) {
        family <- replicate(nresp, family, simplify = FALSE)
    }
    else if (length(family) != nresp) {
        stop2("If 'family' is a list, it has to be of the same ", 
            "length as the number of response variables.")
    }
    if (!is(autocor, "list")) {
        autocor <- replicate(nresp, autocor, simplify = FALSE)
    }
    else if (length(autocor) != nresp) {
        stop2("If 'autocor' is a list, it has to be of the same ", 
            "length as the number of response variables.")
    }
    for (i in seq_len(nresp)) {
        formula$forms[[i]] <- validate_formula(formula$forms[[i]], 
            family = family[[i]], autocor = autocor[[i]], ...)
    }
    if (length(formula$forms) < 2L) {
        stop2("Multivariate models require at least two responses.")
    }
    allow_rescor <- allow_rescor(formula)
    if (is.null(formula$rescor)) {
        miforms <- ulapply(formula$forms, function(f) terms_ad(f$formula, 
            f$family, FALSE)[["mi"]])
        formula$rescor <- allow_rescor && !length(miforms)
        message("Setting 'rescor' to ", formula$rescor, " by default for this model")
        if (formula$rescor) {
            warning2("In the future, 'rescor' will be set to FALSE by default for ", 
                "all models. It is thus recommended to explicitely set ", 
                "'rescor' via 'set_rescor' instead of using the default.")
        }
    }
    formula$rescor <- as_one_logical(formula$rescor)
    if (formula$rescor) {
        if (!allow_rescor) {
            stop2("Currently, estimating 'rescor' is only possible ", 
                "in multivariate gaussian or student models.")
        }
    }
    formula$mecor <- default_mecor(formula$mecor)
    for (i in seq_along(formula$forms)) {
        formula$forms[[i]]$mecor <- formula$mecor
    }
    if (!is.null(cov_ranef)) {
        formula$cov_ranef <- validate_cov_ranef(cov_ranef)
    }
    formula
}
