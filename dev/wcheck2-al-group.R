# Lane wcheck2: does the inner-mode failure of dev/wcheck2-al-decomp.R
# reach an ordinary asym_laplace fit, with an effect per group of 4 rows
# and no weights? The 72-of-80 right-censored data of dev/wcheck-cases.R
# part cens (seed 20261021), yac | cens(ca) ~ x + (1 | g4), and the same
# with every row uncensored. At the fit's own outer parameters: TMB's
# Laplace objective against the Laplace approximation with each group's
# inner mode found by a 1-D maximization, and the exact marginal.
#   Rscript dev/wcheck2-al-group.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$ya <- 1 + 0.5 * d$x + rnorm(n)
o <- order(d$ya, decreasing = TRUE)
d$ca <- "none"
d$ca[o[seq_len(72)]] <- "right"
hi <- d$ya[o[72]]
d$yac <- pmin(d$ya, hi)
d$none <- "none"

row_lp <- function(mu, yi, ci, s, p) {
  u <- (yi - mu) / s
  if (ci == 0) return(log(p * (1 - p) / s) - u * (p - (u < 0)))
  if (u > 0) return(log(1 - p) - p * u)
  log1p(-p * exp((1 - p) * u))
}
row_curv <- function(mu, yi, ci, s, p) {
  u <- (yi - mu) / s
  if (ci == 0 || u > 0) return(0)
  e <- p * exp((1 - p) * u)
  (1 - p)^2 * e / ((1 - e)^2 * s^2)
}
grp <- as.integer(d$g4)
check <- function(lab, fo, cen, yv) {
  fit <- q(frm(fo, family = asym_laplace(), data = d))
  ob <- q(frm(fo, family = asym_laplace(), data = d, dry_run = "objective"))
  obj <- if (!is.null(ob$obj)) ob$obj else ob
  par <- fit$opt$par
  tmb <- -obj$fn(par)
  bh <- obj$env$last.par[obj$env$random]
  nm <- names(par)
  beta <- par[nm == "beta"]
  s <- exp(par[nm == "betad"][1L])
  p <- plogis(par[nm == "betad"][2L])
  sd <- exp(par[nm == "theta"][1L])
  eta <- beta[1L] + beta[2L] * d$x

  f_g <- function(b, k) {
    i <- which(grp == k)
    sum(vapply(i, function(r) row_lp(eta[r] + b, yv[r], cen[r], s, p), 0)) +
      dnorm(b, 0, sd, log = TRUE)
  }
  lap_g <- function(b, k) {
    i <- which(grp == k)
    h <- sum(vapply(i, function(r) row_curv(eta[r] + b, yv[r], cen[r], s, p),
                    0)) + 1 / sd^2
    f_g(b, k) + 0.5 * log(2 * pi) - 0.5 * log(h)
  }
  G <- max(grp)
  tmb_rows <- sum(vapply(seq_len(G), function(k) lap_g(bh[k], k), 0))
  modes <- vapply(seq_len(G), function(k) {
    optimize(function(b) f_g(b, k), c(-50 * sd - 50, 50 * sd + 50),
             maximum = TRUE, tol = 1e-12)$maximum
  }, 0)
  lap <- sum(vapply(seq_len(G), function(k) lap_g(modes[k], k), 0))
  exact <- sum(vapply(seq_len(G), function(k) {
    fm <- f_g(modes[k], k)
    v <- tryCatch(integrate(function(b) exp(vapply(b, f_g, 0, k = k) - fm),
                            modes[k] - 30 * sd, modes[k] + 30 * sd,
                            rel.tol = 1e-10, subdivisions = 5000L)$value,
                  error = function(e) NA_real_)
    fm + log(v)
  }, 0))
  gap <- max(abs(bh - modes))
  cat(sprintf(paste0("%s\n  fit: logLik %.4f conv %d, log sigma %.4f, ",
                     "quantile %.4f, log sd %.4f\n  at the fit: TMB %.4f ",
                     "(rebuilt from TMB's modes %.4f) | exact inner modes ",
                     "%.4f | exact marginal %.4f | max |TMB mode - exact ",
                     "mode| %.3g\n"),
              lab, as.numeric(logLik(fit)), fit$opt$convergence, log(s), p,
              log(sd), tmb, tmb_rows, lap, exact, gap))
}
check("72 of 80 right censored, (1 | g4)",
      bf(yac | cens(ca) ~ x + (1 | g4)), ifelse(d$ca == "right", 1, 0), d$yac)
check("no row censored, (1 | g4)", bf(ya ~ x + (1 | g4)), rep(0, n), d$ya)
cat("DONE\n")
