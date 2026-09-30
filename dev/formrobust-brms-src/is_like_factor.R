is_like_factor <- 
function (x) 
{
    is.factor(x) || is.character(x) || is.logical(x)
}
