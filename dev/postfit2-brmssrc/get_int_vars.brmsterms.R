get_int_vars.brmsterms <- 
function (x, ...) 
{
    adterms <- c("trials", "thres", "vint")
    advars <- ulapply(rmNULL(x$adforms[adterms]), all_vars)
    unique(c(advars, get_sp_vars(x, "mo")))
}
