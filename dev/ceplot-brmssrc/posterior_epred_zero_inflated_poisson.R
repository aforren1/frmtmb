posterior_epred_zero_inflated_poisson <- 
function (prep) 
{
    with(prep$dpars, mu * (1 - zi))
}
