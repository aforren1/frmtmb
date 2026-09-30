conditional_smooths.brmsfit <- 
function (x, smooths = NULL, int_conditions = NULL, prob = 0.95, 
    spaghetti = FALSE, surface = TRUE, resolution = 100, too_far = 0, 
    ndraws = NULL, draw_ids = NULL, nsamples = NULL, subset = NULL, 
    probs = NULL, ...) 
{
    probs <- validate_ci_bounds(prob, probs = probs)
    spaghetti <- as_one_logical(spaghetti)
    surface <- as_one_logical(surface)
    draw_ids <- use_alias(draw_ids, subset)
    ndraws <- use_alias(ndraws, nsamples)
    contains_draws(x)
    x <- restructure(x)
    x <- exclude_terms(x, incl_autocor = FALSE)
    smooths <- rm_wsp(as.character(smooths))
    conditions <- prepare_conditions(x)
    draw_ids <- validate_draw_ids(x, draw_ids, ndraws)
    bterms <- brmsterms(exclude_terms(x$formula, smooths_only = TRUE))
    out <- conditional_smooths(bterms, fit = x, smooths = smooths, 
        conditions = conditions, int_conditions = int_conditions, 
        too_far = too_far, resolution = resolution, probs = probs, 
        spaghetti = spaghetti, surface = surface, draw_ids = draw_ids)
    if (!length(out)) {
        stop2("No valid smooth terms found in the model.")
    }
    structure(out, class = "brms_conditional_effects", smooths_only = TRUE)
}
