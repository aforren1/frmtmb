posterior_epred_wiener <- 
function (prep) 
{
    with(prep$dpars, ndt - bias/mu + bs/mu * (exp(-2 * mu * bias) - 
        1)/(exp(-2 * mu * bs) - 1))
}
