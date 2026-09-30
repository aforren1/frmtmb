prepare_predictions.brmsfit <- 
function (x, newdata = NULL, re_formula = NULL, allow_new_levels = FALSE, 
    sample_new_levels = "uncertainty", incl_autocor = TRUE, oos = NULL, 
    resp = NULL, ndraws = NULL, draw_ids = NULL, nsamples = NULL, 
    subset = NULL, nug = NULL, smooths_only = FALSE, offset = TRUE, 
    newdata2 = NULL, new_objects = NULL, point_estimate = NULL, 
    ndraws_point_estimate = 1, ...) 
{
    x <- restructure(x)
    options(.brmsfit_version = x$version$brms)
    on.exit(options(.brmsfit_version = NULL))
    snl_options <- c("uncertainty", "gaussian", "old_levels")
    sample_new_levels <- match.arg(sample_new_levels, snl_options)
    ndraws <- use_alias(ndraws, nsamples)
    draw_ids <- use_alias(draw_ids, subset)
    warn_brmsfit_multiple(x, newdata = newdata)
    newdata2 <- use_alias(newdata2, new_objects)
    x <- exclude_terms(x, incl_autocor = incl_autocor, offset = offset, 
        smooths_only = smooths_only)
    resp <- validate_resp(resp, x)
    draw_ids <- validate_draw_ids(x, draw_ids, ndraws)
    draws <- as_draws_matrix(x)
    draws <- suppressMessages(subset_draws(draws, draw = draw_ids, 
        unique = FALSE))
    draws <- point_draws(draws, point_estimate, ndraws_point_estimate)
    sdata <- standata(x, newdata = newdata, re_formula = re_formula, 
        newdata2 = newdata2, resp = resp, allow_new_levels = allow_new_levels, 
        internal = TRUE, ...)
    new_formula <- update_re_terms(x$formula, re_formula)
    bframe <- brmsframe(new_formula, data = x$data)
    prep_re <- prepare_predictions_re_global(bframe = bframe, 
        draws = draws, sdata = sdata, resp = resp, old_reframe = x$ranef, 
        sample_new_levels = sample_new_levels, )
    prepare_predictions(bframe, draws = draws, sdata = sdata, 
        prep_re = prep_re, resp = resp, sample_new_levels = sample_new_levels, 
        nug = nug, new = !is.null(newdata), oos = oos, stanvars = x$stanvars)
}
