plot.brmsMarginalEffects <- 
function (x, ...) 
{
    class(x) <- "brms_conditional_effects"
    plot(x, ...)
}
