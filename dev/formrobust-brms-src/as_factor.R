as_factor <- 
function (x, levels = NULL) 
{
    if (is.null(levels)) {
        out <- as.factor(x)
    }
    else {
        out <- factor(x, levels = levels)
    }
    out
}
