posterior_smooths.brmsfit <- 
function (object, smooth, newdata = NULL, resp = NULL, dpar = NULL, 
    nlpar = NULL, ndraws = NULL, draw_ids = NULL, ...) 
{
    resp <- validate_resp(resp, object, multiple = FALSE)
    bterms <- brmsterms(exclude_terms(object$formula, smooths_only = TRUE))
    if (!is.null(resp)) {
        stopifnot(is.mvbrmsterms(bterms))
        bterms <- bterms$terms[[resp]]
    }
    if (!is.null(nlpar)) {
        if (length(dpar)) {
            stop2("Cannot use 'dpar' and 'nlpar' at the same time.")
        }
        nlpar <- as_one_character(nlpar)
        nlpars <- names(bterms$nlpars)
        if (!nlpar %in% nlpars) {
            stop2("Invalid argument 'nlpar'. Valid non-linear ", 
                "parameters are: ", collapse_comma(nlpars))
        }
        bterms <- bterms$nlpars[[nlpar]]
    }
    else {
        dpar <- dpar %||% "mu"
        dpar <- as_one_character(dpar)
        dpars <- names(bterms$dpars)
        if (!dpar %in% dpars) {
            stop2("Invalid argument 'dpar'. Valid distributional ", 
                "parameters are: ", collapse_comma(dpars))
        }
        bterms <- bterms$dpars[[dpar]]
    }
    posterior_smooths(bterms, fit = object, smooth = smooth, 
        newdata = newdata, ndraws = ndraws, draw_ids = draw_ids, 
        ...)
}
