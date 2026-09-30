posterior_epred_poisson <- 
function (prep) 
{
    multiply_dpar_rate_denom(prep$dpars$mu, prep)
}
