posterior_epred_zero_one_inflated_beta <- 
function (prep) 
{
    with(prep$dpars, zoi * coi + mu * (1 - zoi))
}
