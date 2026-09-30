posterior_epred_zero_inflated_negbinomial <- 
function (prep) 
{
    with(prep$dpars, mu * (1 - zi))
}
