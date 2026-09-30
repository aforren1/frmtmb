posterior_epred.brmsfit <- 
function (object, newdata = NULL, re_formula = NULL, re.form = NULL, 
    resp = NULL, dpar = NULL, nlpar = NULL, ndraws = NULL, draw_ids = NULL, 
    sort = FALSE, ...) 
{
    cl <- match.call()
    if ("re.form" %in% names(cl) && !missing(re.form)) {
        re_formula <- re.form
    }
    contains_draws(object)
    object <- restructure(object)
    prep <- prepare_predictions(object, newdata = newdata, re_formula = re_formula, 
        resp = resp, ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, 
        ...)
    posterior_epred(prep, dpar = dpar, nlpar = nlpar, sort = sort, 
        scale = "response", summary = FALSE)
}
