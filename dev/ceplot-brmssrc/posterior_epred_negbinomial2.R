posterior_epred_negbinomial2 <- 
function (prep) 
{
    multiply_dpar_rate_denom(prep$dpars$mu, prep)
}
