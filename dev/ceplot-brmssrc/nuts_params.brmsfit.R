nuts_params.brmsfit <- 
function (object, pars = NULL, ...) 
{
    contains_draws(object)
    bayesplot::nuts_params(object$fit, pars = pars, ...)
}
