posterior_epred_hurdle_gamma <- 
function (prep) 
{
    with(prep$dpars, mu * (1 - hu))
}
