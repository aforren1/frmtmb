posterior_epred_negbinomial <- 
function (prep) 
{
    multiply_dpar_rate_denom(prep$dpars$mu, prep)
}
