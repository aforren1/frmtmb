posterior_epred_zero_inflated_beta_binomial <- 
function (prep) 
{
    trials <- data2draws(prep$data$trials, dim_mu(prep))
    prep$dpars$mu * trials * (1 - prep$dpars$zi)
}
