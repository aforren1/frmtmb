cor_arma <- 
function (formula = ~1, p = 0, q = 0, r = 0, cov = FALSE) 
{
    formula <- as.formula(formula)
    p <- as_one_numeric(p)
    q <- as_one_numeric(q)
    cov <- as_one_logical(cov)
    if ("r" %in% names(match.call())) {
        warning2("The ARR structure is no longer supported and ignored.")
    }
    if (!(p >= 0 && p == round(p))) {
        stop2("Autoregressive order must be a non-negative integer.")
    }
    if (!(q >= 0 && q == round(q))) {
        stop2("Moving-average order must be a non-negative integer.")
    }
    if (!sum(p, q)) {
        stop2("At least one of 'p' and 'q' should be greater zero.")
    }
    if (cov && (p > 1 || q > 1)) {
        stop2("Covariance formulation of ARMA structures is ", 
            "only possible for effects of maximal order one.")
    }
    x <- nlist(formula, p, q, cov)
    class(x) <- c("cor_arma", "cor_brms")
    x
}
