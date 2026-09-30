posterior_epred_binomial <- 
function (prep) 
{
    trials <- data2draws(prep$data$trials, dim_mu(prep))
    prep$dpars$mu * trials
}
