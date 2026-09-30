parnames <- 
function (x, ...) 
{
    warning2("'parnames' is deprecated. Please use 'variables' instead.")
    UseMethod("parnames")
}
