posterior_epred_asym_laplace <- 
function (prep) 
{
    with(prep$dpars, mu + sigma * (1 - 2 * quantile)/(quantile * 
        (1 - quantile)))
}
