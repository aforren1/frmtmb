posterior_epred_hurdle_poisson <- 
function (prep) 
{
    with(prep$dpars, mu/(1 - exp(-mu)) * (1 - hu))
}
