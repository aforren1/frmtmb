posterior_epred_trunc_weibull <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < 0, 0, lb)
    prep$dpars$a <- 1 + 1/prep$dpars$shape
    prep$dpars$scale <- with(prep$dpars, mu/gamma(a))
    m1 <- with(prep$dpars, scale * (incgamma(a, (ub/scale)^shape) - 
        incgamma(a, (lb/scale)^shape)))
    with(prep$dpars, m1/(pweibull(ub, shape, scale = scale) - 
        pweibull(lb, shape, scale = scale)))
}
