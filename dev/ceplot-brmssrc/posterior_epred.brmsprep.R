posterior_epred.brmsprep <- 
function (object, dpar, nlpar, sort, scale = "response", incl_thres = NULL, 
    summary = FALSE, robust = FALSE, probs = c(0.025, 0.975), 
    ...) 
{
    summary <- as_one_logical(summary)
    dpars <- names(object$dpars)
    nlpars <- names(object$nlpars)
    if (length(dpar)) {
        dpar <- as_one_character(dpar)
        if (!dpar %in% dpars) {
            stop2("Invalid argument 'dpar'. Valid distributional ", 
                "parameters are: ", collapse_comma(dpars))
        }
        if (length(nlpar)) {
            stop2("Cannot use 'dpar' and 'nlpar' at the same time.")
        }
        predicted <- is.bprepl(object$dpars[[dpar]]) || is.bprepnl(object$dpars[[dpar]])
        if (predicted) {
            if (scale == "linear") {
                object$dpars[[dpar]]$family$link <- "identity"
            }
            if (is_ordinal(object$family)) {
                object$dpars[[dpar]]$cs <- NULL
                object$family <- object$dpars[[dpar]]$family <- .dpar_family(link = object$dpars[[dpar]]$family$link)
            }
            if (dpar_class(dpar) == "theta" && scale == "response") {
                ap_id <- as.numeric(dpar_id(dpar))
                out <- get_theta(object)[, , ap_id, drop = FALSE]
                dim(out) <- dim(out)[c(1, 2)]
            }
            else {
                out <- get_dpar(object, dpar = dpar, inv_link = TRUE)
            }
        }
        else {
            out <- object$dpars[[dpar]]
            out <- matrix(out, nrow = object$ndraws, ncol = object$nobs)
        }
    }
    else if (length(nlpar)) {
        nlpar <- as_one_character(nlpar)
        if (!nlpar %in% nlpars) {
            stop2("Invalid argument 'nlpar'. Valid non-linear ", 
                "parameters are: ", collapse_comma(nlpars))
        }
        out <- get_nlpar(object, nlpar = nlpar)
    }
    else {
        incl_thres <- as_one_logical(incl_thres %||% FALSE)
        incl_thres <- incl_thres && is_ordinal(object$family) && 
            scale == "linear"
        if (incl_thres) {
            if (is.mixfamily(object$family)) {
                stop2("'incl_thres' is not supported for mixture models.")
            }
            object$family$link <- "identity"
        }
        if (scale == "response" || incl_thres) {
            for (nlp in nlpars) {
                object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
            }
            for (dp in dpars) {
                object$dpars[[dp]] <- get_dpar(object, dpar = dp)
            }
            if (is_trunc(object)) {
                out <- posterior_epred_trunc(object)
            }
            else {
                posterior_epred_fun <- paste0("posterior_epred_", 
                  object$family$family)
                posterior_epred_fun <- get(posterior_epred_fun, 
                  asNamespace("brms"))
                out <- posterior_epred_fun(object)
            }
        }
        else {
            if (conv_cats_dpars(object$family)) {
                out <- dpars[grepl("^mu", dpars)]
            }
            else {
                out <- dpars[dpar_class(dpars) %in% "mu"]
            }
            if (length(out) == 1) {
                out <- get_dpar(object, dpar = out, inv_link = FALSE)
            }
            else {
                out <- lapply(out, get_dpar, prep = object, inv_link = FALSE)
                out <- abind::abind(out, along = 3)
            }
        }
    }
    if (is.null(dim(out))) {
        out <- as.matrix(out)
    }
    colnames(out) <- NULL
    out <- reorder_obs(out, object$old_order, sort = sort)
    if (scale == "response" && is_polytomous(object$family) && 
        length(dim(out)) == 3 && dim(out)[3] == length(object$cats)) {
        dimnames(out)[[3]] <- object$cats
    }
    if (summary) {
        out <- posterior_summary(out, probs = probs, robust = robust)
        if (is_polytomous(object$family) && length(dim(out)) == 
            3) {
            if (scale == "linear") {
                dimnames(out)[[3]] <- paste0("eta", seq_dim(out, 
                  3))
            }
            else {
                dimnames(out)[[3]] <- paste0("P(Y = ", dimnames(out)[[3]], 
                  ")")
            }
        }
    }
    out
}
