# Review of lane alapse: why ald_vertex() fails on tied data. Calls it
# on the ties cases of alapse-review-rq.R at rq's own solution, and
# counts zero residuals there (a degenerate vertex has more than k).
#   Rscript dev/reviews/alapse-review-ties.R lane
lib_arg <- "lane"
source("dev/reviews/alapse-review-common.R")
for (s in 1:3) {
  set.seed(20 + s)
  n <- 60
  d <- data.frame(x = sample(0:2, n, TRUE))
  d$y <- rpois(n, 3 + d$x)
  X <- cbind(1, d$x)
  rr <- quantreg::rq.fit(X, d$y, tau = 0.5, method = "br")
  b <- rr$coefficients
  res <- d$y - drop(X %*% b)
  cat(sprintf("s%d: rq b = %s, zero residuals %d (k = 2)\n", s,
              paste(format(b), collapse = " "), sum(abs(res) < 1e-12)))
  for (b0 in list(b, b + c(0.01, -0.01), c(mean(d$y), 0))) {
    v <- ns$ald_vertex(X, d$y, rep(1, n), rep(0.5, n), b0)
    cat("  from", format(b0, digits = 4), ":",
        if (is.null(v)) "NULL" else sprintf("b %s margin %.3g iter %d",
                                            paste(format(v$b), collapse = " "),
                                            v$margin, v$iter), "\n")
    v2 <- ns$ald_vertex(X, d$y, rep(1, n), rep(0.5, n), b0, iter = 5000L)
    cat("    iter 5000:", if (is.null(v2)) "NULL" else
      sprintf("b %s margin %.3g iter %d", paste(format(v2$b), collapse = " "),
              v2$margin, v2$iter), "\n")
  }
  # the same data with ties broken by a jitter of 1e-9
  set.seed(1)
  yj <- d$y + runif(n, -1e-9, 1e-9)
  v3 <- ns$ald_vertex(X, yj, rep(1, n), rep(0.5, n), b)
  cat("  jittered:", if (is.null(v3)) "NULL" else
    sprintf("b %s margin %.3g iter %d", paste(format(v3$b), collapse = " "),
            v3$margin, v3$iter), "\n")
}
cat("DONE\n")
