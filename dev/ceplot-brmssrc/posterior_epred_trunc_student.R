posterior_epred_trunc_student <- 
function (prep, lb, ub) 
{
    zlb <- with(prep$dpars, (lb - mu)/sigma)
    zub <- with(prep$dpars, (ub - mu)/sigma)
    nu <- prep$dpars$nu
    G1 <- gamma((nu - 1)/2) * nu^(nu/2)/(2 * (pt(zub, df = nu) - 
        pt(zlb, df = nu)) * gamma(nu/2) * gamma(0.5))
    A <- (nu + zlb^2)^(-(nu - 1)/2)
    B <- (nu + zub^2)^(-(nu - 1)/2)
    trunc_zmean <- G1 * (A - B)
    prep$dpars$mu + trunc_zmean * prep$dpars$sigma
}
