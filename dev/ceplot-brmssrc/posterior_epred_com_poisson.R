posterior_epred_com_poisson <- 
function (prep) 
{
    mean_com_poisson(prep$dpars$mu, prep$dpars$shape)
}
