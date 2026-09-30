nwarmup <- 
function (x) 
{
    if (!is.stanfit(x$fit)) 
        return(0)
    x$fit@sim$warmup2[1] %||% 0
}
