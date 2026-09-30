posterior_epred_shifted_lognormal <- 
function (prep) 
{
    with(prep$dpars, exp(mu + sigma^2/2) + ndt)
}
