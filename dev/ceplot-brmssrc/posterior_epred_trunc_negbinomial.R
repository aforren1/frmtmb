posterior_epred_trunc_negbinomial <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < -1, -1, lb)
    mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
    max_value <- 3 * max(mu)
    ub <- ifelse(ub > max_value, max_value, ub)
    shape <- multiply_dpar_rate_denom(prep$dpars$shape, prep)
    args <- list(mu = mu, size = shape)
    posterior_epred_trunc_discrete(dist = "nbinom", args = args, 
        lb = lb, ub = ub)
}
