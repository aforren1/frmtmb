posterior_epred_gen_extreme_value <- 
function (prep) 
{
    with(prep$dpars, mu + sigma * (gamma(1 - xi) - 1)/xi)
}
