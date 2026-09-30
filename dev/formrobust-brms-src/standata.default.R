standata.default <- 
function (object, data, family = gaussian(), prior = NULL, autocor = NULL, 
    data2 = NULL, cov_ranef = NULL, sample_prior = "no", stanvars = NULL, 
    threads = getOption("brms.threads", NULL), knots = NULL, 
    drop_unused_levels = TRUE, ...) 
{
    object <- validate_formula(object, data = data, family = family, 
        autocor = autocor, cov_ranef = cov_ranef)
    bterms <- brmsterms(object)
    data2 <- validate_data2(data2, bterms = bterms, get_data2_autocor(object), 
        get_data2_cov_ranef(object))
    data <- validate_data(data, bterms = bterms, knots = knots, 
        data2 = data2, drop_unused_levels = drop_unused_levels)
    bframe <- brmsframe(bterms, data)
    prior <- .validate_prior(prior, bframe = bframe, sample_prior = sample_prior)
    stanvars <- validate_stanvars(stanvars)
    threads <- validate_threads(threads)
    .standata(bframe, data = data, prior = prior, data2 = data2, 
        stanvars = stanvars, threads = threads, ...)
}
