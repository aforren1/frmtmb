posterior_epred_lagsar <- 
function (prep) 
{
    stopifnot(!is.null(prep$ac$lagsar))
    I <- diag(prep$nobs)
    .posterior_epred <- function(s) {
        IB <- I - with(prep$ac, lagsar[s, ] * Msar)
        as.numeric(solve(IB, prep$dpars$mu[s, ]))
    }
    out <- rblapply(seq_len(prep$ndraws), .posterior_epred)
    rownames(out) <- NULL
    out
}
