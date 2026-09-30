fitted.brmsfit <- 
function (object, newdata = NULL, re_formula = NULL, scale = c("response", 
    "linear"), resp = NULL, dpar = NULL, nlpar = NULL, ndraws = NULL, 
    draw_ids = NULL, sort = FALSE, summary = TRUE, robust = FALSE, 
    probs = c(0.025, 0.975), ...) 
{
    scale <- match.arg(scale)
    summary <- as_one_logical(summary)
    contains_draws(object)
    object <- restructure(object)
    prep <- prepare_predictions(object, newdata = newdata, re_formula = re_formula, 
        resp = resp, ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, 
        ...)
    posterior_epred(prep, dpar = dpar, nlpar = nlpar, sort = sort, 
        scale = scale, summary = summary, robust = robust, probs = probs)
}
