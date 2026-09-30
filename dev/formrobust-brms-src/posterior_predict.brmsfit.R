posterior_predict.brmsfit <- 
function (object, newdata = NULL, re_formula = NULL, re.form = NULL, 
    transform = NULL, resp = NULL, negative_rt = FALSE, ndraws = NULL, 
    draw_ids = NULL, sort = FALSE, ntrys = 5, cores = NULL, ...) 
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
    posterior_predict(prep, transform = transform, sort = sort, 
        ntrys = ntrys, negative_rt = negative_rt, cores = cores, 
        summary = FALSE)
}
