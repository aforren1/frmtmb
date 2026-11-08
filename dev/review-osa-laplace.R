# Review lane osa: is the gap to the exact one-step predictive the Laplace
# approximation, or the route? Poisson, (1 | g), the seed-201 design of
# dev/review-osa-exact.R. The reviewer's own Laplace approximation of each
# of TMB's three integrals (density, lower tail, upper tail), in R, from
# dpois()/ppois(); then a quadrature = TRUE fit of the same data.
#
#   Rscript dev/review-osa-laplace.R [seed]
args <- commandArgs(TRUE)
seed <- if (length(args)) as.integer(args[1]) else 201L
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-laplace seed", seed, "\n")
G <- 20L
K <- 10L
n <- G * K
U <- local({ set.seed(123); runif(n) })
set.seed(seed)
d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                x = rnorm(n), k = sample(1:12, n, replace = TRUE))
e <- rnorm(G, 0, 1.0)[d$g]
d$y <- rpois(n, exp(0.2 + 0.5 * d$x + e))
fit <- suppressWarnings(frm(y ~ x + (1 | g), family = poisson(), data = d))
r <- residuals(fit, type = "osa")[, "Estimate"]
eta <- fit$estimates$beta[[1]] + fit$estimates$beta[[2]] * d$x
sdu <- exp(fit$estimates$theta[[1]])

lap <- function(h) {
  o <- optimize(function(u) -h(u), c(-10, 10) * sdu)
  u0 <- o$minimum
  e <- 1e-4 * max(1, abs(u0))
  h2 <- (h(u0 + e) - 2 * h(u0) + h(u0 - e)) / e^2
  h(u0) + 0.5 * log(2 * pi) - 0.5 * log(-h2)
}
quad <- function(h) {
  u <- seq(-8, 8, length.out = 1601) * sdu
  v <- h(u)
  m <- max(v)
  m + log(sum(exp(v - m)) * (u[2] - u[1]))
}
own <- function(integ) {
  out <- numeric(n)
  for (gi in split(seq_len(n), d$g)) {
    for (k in seq_along(gi)) {
      i <- gi[k]
      prev <- gi[seq_len(k - 1)]
      base <- function(u) {
        s <- dnorm(u, 0, sdu, log = TRUE)
        for (j in prev) s <- s + dpois(d$y[j], exp(eta[j] + u), log = TRUE)
        s
      }
      lp <- integ(function(u) base(u) + dpois(d$y[i], exp(eta[i] + u), log = TRUE))
      lF <- integ(function(u) base(u) + ppois(d$y[i], exp(eta[i] + u), log.p = TRUE))
      lS <- integ(function(u) base(u) + ppois(d$y[i], exp(eta[i] + u),
                                               lower.tail = FALSE, log.p = TRUE))
      Fx <- 1 / (1 + exp(lS - lF))
      px <- exp(lp - lF) * Fx
      out[i] <- qnorm(Fx - U[i] * px)
    }
  }
  out
}
rl <- own(lap)
rq <- own(quad)
first <- d$t == 1
cat(sprintf("lane vs own Laplace:  max %.3g median %.3g\n",
            max(abs(r - rl)), median(abs(r - rl))))
cat(sprintf("lane vs exact:        max %.3g median %.3g; first rows of a group max %.3g, later rows max %.3g\n",
            max(abs(r - rq)), median(abs(r - rq)), max(abs(r - rq)[first]),
            max(abs(r - rq)[!first])))
cat(sprintf("own Laplace vs exact: max %.3g median %.3g\n",
            max(abs(rl - rq)), median(abs(rl - rq))))
fq <- tryCatch(suppressWarnings(frm(y ~ x + (1 | g), family = poisson(),
                                    data = d, quadrature = TRUE)),
               error = function(e) e)
if (inherits(fq, "error")) {
  cat("quadrature = TRUE fit:", conditionMessage(fq), "\n")
} else {
  rqf <- tryCatch(residuals(fq, type = "osa")[, "Estimate"],
                  error = function(e) e)
  if (inherits(rqf, "error")) {
    cat("quadrature = TRUE osa:", conditionMessage(rqf), "\n")
  } else {
    cat(sprintf(paste("quadrature = TRUE fit: theta %.5f vs Laplace %.5f;",
                      "osa vs Laplace-fit osa max %.3g; the osa is still",
                      "the Laplace one-step predictive at the quadrature",
                      "estimates\n"),
                fq$estimates$theta[[1]], fit$estimates$theta[[1]],
                max(abs(rqf - r))))
  }
}
cat("END review-osa-laplace\n")
