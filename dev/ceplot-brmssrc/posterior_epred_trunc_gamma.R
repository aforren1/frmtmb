posterior_epred_trunc_gamma <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < 0, 0, lb)
    prep$dpars$scale <- prep$dpars$mu/prep$dpars$shape
    m1 <- with(prep$dpars, scale/gamma(shape) * (incgamma(1 + 
        shape, ub/scale) - incgamma(1 + shape, lb/scale)))
    with(prep$dpars, m1/(pgamma(ub, shape, scale = scale) - pgamma(lb, 
        shape, scale = scale)))
}
