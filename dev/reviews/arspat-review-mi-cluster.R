# Review of lane arspat: vcov_cluster() on mi() fits. Docs list mi() as
# refused; check a plain mi() fit (before and after) and the ar() + mi()
# fit against an independent sandwich of the dense margin.
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(numDeriv))
set.seed(17)
G <- 30; K <- 6
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n) + rep(rnorm(G), each = K)
e <- as.vector(vapply(seq_len(G), function(i)
  as.vector(arima.sim(list(ar = 0.5), K, sd = 0.7)), numeric(K)))
d$y <- 1 + 0.6 * d$x + e
miss <- c(2, 3, 13, 40, 41, 42, 90, 120, 179)
d$y[miss] <- NA
chk("plain y | mi() ~ x, vcov_cluster(~ g)", {
  f <- frm(y | mi() ~ x, data = d)
  vcov_cluster(f, ~ g, type = "CR0")
})
chk("ar + mi, vcov_cluster against the dense margin", {
  f <- frm(y | mi() ~ x + ar(t, g, cov = TRUE), data = d)
  V <- vcov_cluster(f, ~ g, type = "CR0")
  ll <- function(th, i) {
    sig <- exp(th[3]); phi <- tanh(th[4])
    i <- i[!is.na(d$y[i])]
    tt <- d$t[i]
    S <- sig^2 / (1 - phi^2) * phi^abs(outer(tt, tt, "-"))
    mvtnorm::dmvnorm(d$y[i], th[1] + th[2] * d$x[i], S, log = TRUE)
  }
  grp <- split(seq_len(n), d$g)
  tot <- function(th) sum(vapply(grp, function(i) ll(th, i), 0))
  o <- optim(c(0, 0, 0, 0), function(th) -tot(th), method = "BFGS",
             control = list(reltol = 1e-14, maxit = 1000))
  H <- hessian(tot, o$par)
  S <- t(vapply(grp, function(i) grad(function(t) ll(t, i), o$par),
                numeric(4)))
  B <- solve(-H)
  Vr <- (B %*% crossprod(S) %*% B)[1:2, 1:2]
  cat("V frmtmb:", format(as.vector(V[1:2, 1:2]), digits = 8), "\n")
  cat("V indep :", format(as.vector(Vr), digits = 8), "\n")
  NULL
})
