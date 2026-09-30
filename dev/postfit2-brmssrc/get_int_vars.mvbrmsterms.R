get_int_vars.mvbrmsterms <- 
function (x, ...) 
{
    unique(ulapply(x$terms, get_int_vars))
}
