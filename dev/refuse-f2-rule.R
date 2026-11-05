# N5: which rule integrates a read value with the fewest sweeps of the objective
# and the same total: integrate() on the infinite range, against one finite range
# through z = w / (1 - w^2). Reference: a core and two tails at rel.tol 1e-6.
# Seed 11503.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
pr <- function(...) cat(sprintf(...), "\n", sep = "")
rows_err <- function(fit, rule) {
  env <- fit$obj$env; par <- env$last.par.best
  mp <- which(names(par) == "miss")
  H <- env$spHess(par, random = TRUE)
  hk <- Matrix::diag(H)[match(mp, env$random)]
  f0 <- env$f(par); count <- 0L
  chart <- character(length(mp))
  for (m in fit$frame$mi_map) chart[m$idx] <- m$chart %||% "identity"
  err <- vapply(seq_along(mp), function(k) {
    u0 <- par[mp[k]]; s <- 1 / sqrt(hk[k])
    h <- function(z) vapply(z, function(zi) {
      count <<- count + 1L
      p <- par; p[mp[k]] <- u0 + s * zi
      v <- exp(f0 - env$f(p)); if (is.finite(v)) v else 0
    }, 0)
    i <- switch(rule,
      inf = stats::integrate(h, -Inf, Inf, subdivisions = 100L, rel.tol = 1e-2, stop.on.error = FALSE)$value,
      inf10 = stats::integrate(h, -Inf, Inf, subdivisions = 100L, rel.tol = 0.1, stop.on.error = FALSE)$value,
      w = stats::integrate(function(w) h(w / (1 - w^2)) * (1 + w^2) / (1 - w^2)^2, -1, 1,
                           subdivisions = 100L, rel.tol = 1e-2, stop.on.error = FALSE)$value,
      ref = {
        lo <- if (chart[k] %in% c("logit", "circle")) (-37 - u0) / s else -Inf
        hi <- if (chart[k] %in% c("logit", "circle")) (37 - u0) / s else Inf
        pa <- function(a, b) stats::integrate(h, a, b, subdivisions = 400L, rel.tol = 1e-6, stop.on.error = FALSE)$value
        pa(max(-6, lo), min(6, hi)) + pa(lo, max(-6, lo)) + pa(min(6, hi), hi)
      })
    0.5 * log(2 * pi) - log(i)
  }, 0)
  c(total = sum(err), sweeps = count / length(mp))
}
rvm <- function(mu, k) {
  vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
}
set.seed(11503)
n <- 200; x <- rnorm(n)
cases <- list(
  exg_weak = list(fam = exgaussian(), y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2, b = 0.05),
  exg_strong = list(fam = exgaussian(), y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2, b = 0.7),
  gamma = list(fam = stats::Gamma(link = "log"), y = rgamma(n, 1.3, 1.3 / exp(0.3 + 0.5 * x)), b = 0.7),
  beta = list(fam = Beta(), y = rbeta(n, plogis(0.3 + 0.7 * x) * 3, (1 - plogis(0.3 + 0.7 * x)) * 3), b = 0.7),
  student = list(fam = student(), y = 0.3 + 0.5 * x + rt(n, 3), b = 0.05),
  vm_mid = list(fam = von_mises(), y = rvm(2 * atan(0.3 + 0.5 * x), 3), b = 0.7),
  vm_end = list(fam = von_mises(), y = rvm(3.05 + 0.05 * x, 3), b = 0.5),
  vm_end20 = list(fam = von_mises(), y = rvm(3.14 + 0.05 * x, 20), b = 0.5))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- data.frame(x = x, y = cs$y); d$z <- 0.5 + cs$b * cs$y + rnorm(n, 0, 0.7)
  d$y[sample(n, 40)] <- NA
  fit <- suppressWarnings(frm(bf(y | mi() ~ x) + cs$fam + bf(z ~ mi(y)) + gaussian(), data = d))
  a <- rows_err(fit, "ref"); b <- rows_err(fit, "inf"); w <- rows_err(fit, "inf10")
  pr("%-10s ref %8.3f (%4.0f sweeps); inf %8.3f (%4.0f); inf at 0.1 %8.3f (%4.0f)", nm, a[1], a[2], b[1], b[2], w[1], w[2])
}
