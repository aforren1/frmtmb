function (x) 
{
    no_int <- no_int(x)
    no_cmc <- no_cmc(x)
    if (is.formula(x) && !is.terms(x)) {
        x <- terms(x)
    }
    if (!is.terms(x)) {
        return(NULL)
    }
    if (no_int || !has_intercept(x) && no_cmc) {
        attr(x, "intercept") <- 1
        attr(x, "int") <- FALSE
    }
    x
}
