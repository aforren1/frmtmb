as_formula <- 
function (x) 
{
    x <- as.formula(x)
    rhs <- rhs(x)[[2]]
    if (isTRUE(is.call(rhs) && rhs[[1]] == "~")) {
        stop2("Nested formulas are not allowed. Did you use '~~' somewhere?")
    }
    x
}
