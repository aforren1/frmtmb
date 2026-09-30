predictive_error.brmsfit <- 
function (object, newdata = NULL, re_formula = NULL, re.form = NULL, 
    method = "posterior_predict", resp = NULL, ndraws = NULL, 
    draw_ids = NULL, sort = FALSE, ...) 
{
    cl <- match.call()
    if ("re.form" %in% names(cl) && !missing(re.form)) {
        re_formula <- re.form
    }
    .predictive_error(object, newdata = newdata, re_formula = re_formula, 
        method = method, type = "ordinary", resp = resp, ndraws = ndraws, 
        draw_ids = draw_ids, sort = sort, ...)
}
