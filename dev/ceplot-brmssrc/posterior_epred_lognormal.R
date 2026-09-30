posterior_epred_lognormal <- 
function (prep) 
{
    with(prep$dpars, exp(mu + sigma^2/2))
}
