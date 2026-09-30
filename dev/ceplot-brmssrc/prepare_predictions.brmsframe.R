prepare_predictions.brmsframe <- 
function (x, draws, sdata, ...) 
{
    ndraws <- nrow(draws)
    nobs <- sdata[[paste0("N", usc(x$resp))]]
    resp <- usc(combine_prefix(x))
    out <- nlist(ndraws, nobs, resp = x$resp)
    out$family <- prepare_family(x)
    out$old_order <- attr(sdata, "old_order")
    if (has_subset(x) && !is.null(out$old_order)) {
        out$old_order <- as.numeric(factor(out$old_order[x$frame$resp$subset]))
    }
    valid_dpars <- valid_dpars(x)
    out$dpars <- named_list(valid_dpars)
    for (dp in valid_dpars) {
        dp_regex <- paste0("^", dp, resp, "$")
        if (is.btl(x$dpars[[dp]]) || is.btnl(x$dpars[[dp]])) {
            out$dpars[[dp]] <- prepare_predictions(x$dpars[[dp]], 
                draws = draws, sdata = sdata, ...)
        }
        else if (any(grepl(dp_regex, colnames(draws)))) {
            out$dpars[[dp]] <- as.vector(prepare_draws(draws, 
                dp_regex, regex = TRUE))
        }
        else if (is.numeric(x$fdpars[[dp]]$value)) {
            out$dpars[[dp]] <- x$fdpars[[dp]]$value
        }
    }
    out$nlpars <- named_list(names(x$nlpars))
    for (nlp in names(x$nlpars)) {
        out$nlpars[[nlp]] <- prepare_predictions(x$nlpars[[nlp]], 
            draws = draws, sdata = sdata, ...)
    }
    if (is.mixfamily(x$family)) {
        families <- family_names(x$family)
        thetas <- paste0("theta", seq_along(families))
        if (any(ulapply(out$dpars[thetas], is.list))) {
            missing_id <- which(ulapply(out$dpars[thetas], is.null))
            out$dpars[[paste0("theta", missing_id)]] <- structure(data2draws(0, 
                c(ndraws, nobs)), predicted = TRUE)
        }
        else {
            out$dpars$theta <- do_call(cbind, out$dpars[thetas])
            out$dpars[thetas] <- NULL
            if (nrow(out$dpars$theta) == 1) {
                dim <- c(nrow(draws), ncol(out$dpars$theta))
                out$dpars$theta <- data2draws(out$dpars$theta, 
                  dim = dim)
            }
        }
    }
    if (is_ordinal(x$family)) {
        if (is.mixfamily(x$family)) {
            mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
            for (mu in mu_pars) {
                out$thres[[mu]] <- prepare_predictions_thres(x$dpars[[mu]], 
                  draws, sdata, ...)
            }
        }
        else {
            out$thres <- prepare_predictions_thres(x$dpars$mu, 
                draws, sdata, ...)
        }
    }
    if (is_logistic_normal(x$family)) {
        out$dpars$lncor <- prepare_draws(draws, "^lncor__", regex = TRUE)
    }
    if (is_cox(x$family)) {
        if (is.mixfamily(x$family)) {
            mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
            for (mu in mu_pars) {
                out$bhaz[[mu]] <- prepare_predictions_bhaz(x$dpars[[mu]], 
                  draws, sdata, ...)
            }
        }
        else {
            out$bhaz <- prepare_predictions_bhaz(x$dpars$mu, 
                draws, sdata, ...)
        }
    }
    out$cats <- get_cats(x)
    out$refcat <- get_refcat(x, int = TRUE)
    out$ac <- prepare_predictions_ac(x$dpars$mu, draws, sdata, 
        nat_cov = TRUE, ...)
    out$data <- prepare_predictions_data(x, sdata = sdata, ...)
    structure(out, class = "brmsprep")
}
