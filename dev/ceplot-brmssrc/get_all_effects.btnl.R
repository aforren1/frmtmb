get_all_effects.btnl <- 
function (x, ...) 
{
    covars <- all_vars(rhs(x$covars))
    out <- as.list(covars)
    if (length(covars) > 1) {
        c(out) <- utils::combn(covars, 2, simplify = FALSE)
    }
    unique(out)
}
