posterior_samples <- 
function (x, pars = NA, ...) 
{
    warning2("Method 'posterior_samples' is deprecated. ", "Please see ?as_draws for recommended alternatives.")
    UseMethod("posterior_samples")
}
