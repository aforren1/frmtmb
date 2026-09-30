recover_data.brmsfit <- 
function (object, data, resp = NULL, dpar = NULL, nlpar = NULL, 
    re_formula = NA, epred = FALSE, ...) 
{
    bterms <- .extract_par_terms(object, resp = resp, dpar = dpar, 
        nlpar = nlpar, re_formula = re_formula, epred = epred)
    data <- rm_attr(object$data, "terms")
    mf <- model.frame(bterms$allvars, data = data)
    trms <- attr(mf, "terms")
    cl <- call("brms")
    if (epred) {
        cl$formula <- bterms$respform
    }
    emmeans::recover_data(cl, trms, "na.omit", data = data, ...)
}
