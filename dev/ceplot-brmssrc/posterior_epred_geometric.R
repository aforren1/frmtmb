posterior_epred_geometric <- 
function (prep) 
{
    multiply_dpar_rate_denom(prep$dpars$mu, prep)
}
