validate_weights_method <- 
function (method) 
{
    method <- as_one_character(method)
    method <- tolower(method)
    if (method == "loo2") {
        warning2("Weight method 'loo2' is deprecated. Use 'stacking' instead.")
        method <- "stacking"
    }
    if (method == "marglik") {
        warning2("Weight method 'marglik' is deprecated. Use 'bma' instead.")
        method <- "bma"
    }
    options <- c("loo", "waic", "kfold", "stacking", "pseudobma", 
        "bma")
    match.arg(method, options)
}
