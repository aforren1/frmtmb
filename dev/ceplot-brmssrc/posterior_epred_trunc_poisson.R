posterior_epred_trunc_poisson <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < -1, -1, lb)
    mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
    max_value <- 3 * max(mu)
    ub <- ifelse(ub > max_value, max_value, ub)
    args <- list(lambda = mu)
    posterior_epred_trunc_discrete(dist = "pois", args = args, 
        lb = lb, ub = ub)
}
