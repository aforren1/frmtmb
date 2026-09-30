posterior_epred_trunc_lognormal <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < 0, 0, lb)
    m1 <- with(prep$dpars, exp(mu + sigma^2/2) * (pnorm((log(ub) - 
        mu)/sigma - sigma) - pnorm((log(lb) - mu)/sigma - sigma)))
    with(prep$dpars, m1/(plnorm(ub, meanlog = mu, sdlog = sigma) - 
        plnorm(lb, meanlog = mu, sdlog = sigma)))
}
