posterior_epred_trunc_binomial <- 
function (prep, lb, ub) 
{
    lb <- ifelse(lb < -1, -1, lb)
    max_value <- max(prep$data$trials)
    ub <- ifelse(ub > max_value, max_value, ub)
    trials <- prep$data$trials
    if (length(trials) > 1) {
        trials <- data2draws(trials, dim_mu(prep))
    }
    args <- list(size = trials, prob = prep$dpars$mu)
    posterior_epred_trunc_discrete(dist = "binom", args = args, 
        lb = lb, ub = ub)
}
