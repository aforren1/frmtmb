all_vars <- 
function (expr, ...) 
{
    if (is.character(expr)) {
        expr <- str2expression(expr)
    }
    all.vars(expr, ...)
}
