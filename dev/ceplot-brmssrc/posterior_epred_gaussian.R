posterior_epred_gaussian <- 
function (prep) 
{
    if (!is.null(prep$ac$lagsar)) {
        prep$dpars$mu <- posterior_epred_lagsar(prep)
    }
    prep$dpars$mu
}
