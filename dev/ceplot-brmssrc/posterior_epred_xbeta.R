posterior_epred_xbeta <- 
function (prep) 
{
    mu <- get_dpar(prep, "mu")
    phi <- get_dpar(prep, "phi")
    nu <- get_dpar(prep, "kappa")
    a <- mu * phi
    b <- (1 - mu) * phi
    d <- (1 + 2 * nu)
    q0 <- nu/d
    q1 <- (1 + nu)/d
    t3 <- pbeta(q1, a, b)
    t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
    t2 <- nu * (t3 - pbeta(q0, a, b))
    1 + t1 - t2 - t3
}
