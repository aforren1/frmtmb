posterior_epred_trunc_gaussian <- 
function (prep, lb, ub) 
{
    zlb <- (lb - prep$dpars$mu)/prep$dpars$sigma
    zub <- (ub - prep$dpars$mu)/prep$dpars$sigma
    trunc_zmean <- (dnorm(zlb) - dnorm(zub))/(pnorm(zub) - pnorm(zlb))
    prep$dpars$mu + trunc_zmean * prep$dpars$sigma
}
