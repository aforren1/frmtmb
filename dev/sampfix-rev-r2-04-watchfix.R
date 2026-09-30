# Reviewer re-check: where draws_laplace_watch()'s time goes, and the
# same check with unlist(use.names = FALSE). 1000 calls, 2 x 20000 cells.
set.seed(5)
dp <- list(mu = rnorm(20000), sigma = rep(0.7, 20000))
mk <- function(use_names) {
  first <- NULL
  function(v) {
    nf <- !is.finite(suppressWarnings(as.numeric(unlist(v, use.names = use_names))))
    if (is.null(first)) first <<- nf else if (!identical(nf, first)) stop("x")
    invisible(NULL)
  }
}
for (u in c(TRUE, FALSE)) {
  w <- mk(u)
  t <- system.time(for (i in 1:1000) w(dp))[["elapsed"]]
  cat(sprintf("use.names = %-5s: %.3f s for 1000 draws\n", u, t))
}
