posterior_epred_discrete_weibull <- 
function (prep) 
{
    mean_discrete_weibull(prep$dpars$mu, prep$dpars$shape)
}
