posterior_epred_hurdle_lognormal <- 
function (prep) 
{
    with(prep$dpars, exp(mu + sigma^2/2) * (1 - hu))
}
