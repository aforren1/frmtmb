posterior_epred_zero_inflated_asym_laplace <- 
function (prep) 
{
    posterior_epred_asym_laplace(prep) * (1 - prep$dpars$zi)
}
