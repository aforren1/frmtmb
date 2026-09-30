posterior_epred_hurdle_negbinomial <- 
function (prep) 
{
    with(prep$dpars, mu/(1 - (shape/(mu + shape))^shape) * (1 - 
        hu))
}
