posterior_summary.default <- 
function (x, probs = c(0.025, 0.975), robust = FALSE, ...) 
{
    if (!length(x)) {
        stop2("No posterior draws supplied.")
    }
    if (robust) {
        coefs <- c("median", "mad", "quantile")
    }
    else {
        coefs <- c("mean", "sd", "quantile")
    }
    .posterior_summary <- function(x) {
        do_call(cbind, lapply(coefs, get_estimate, draws = x, 
            probs = probs, na.rm = TRUE))
    }
    if (length(dim(x)) <= 2) {
        x <- as.matrix(x)
    }
    else {
        x <- as.array(x)
    }
    if (length(dim(x)) == 2) {
        out <- .posterior_summary(x)
        rownames(out) <- colnames(x)
    }
    else if (length(dim(x)) == 3) {
        out <- lapply(array2list(x), .posterior_summary)
        out <- abind(out, along = 3)
        dnx <- dimnames(x)
        dimnames(out) <- list(dnx[[2]], dimnames(out)[[2]], dnx[[3]])
    }
    else {
        stop("'x' must be of dimension 2 or 3.")
    }
    colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 
        100))
    out
}
