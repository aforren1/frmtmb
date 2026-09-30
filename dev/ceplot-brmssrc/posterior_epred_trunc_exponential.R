posterior_epred_trunc_exponential <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < 0, 0, lb)
    inv_mu <- 1/prep$dpars$mu
    m1 <- with(prep$dpars, mu * (incgamma(2, ub/mu) - incgamma(2, 
        lb/mu)))
    m1/(pexp(ub, rate = inv_mu) - pexp(lb, rate = inv_mu))
}
