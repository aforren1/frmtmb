# What the guard and the domain check cost: nl_check_body_read() and
# nl_domain_start() against a whole frm() on three nl models. Each is timed
# in blocks past 1.2 s, minimum of 5 rounds, with the control of a block of
# the same fits on the same process. Seed 8101.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
block <- function(f) {
  best <- Inf
  for (r in 1:5) {
    k <- 1L
    repeat {
      t0 <- proc.time()[["elapsed"]]
      for (i in seq_len(k)) f()
      dt <- proc.time()[["elapsed"]] - t0
      if (dt > 1.2) break
      k <- k * 2L
    }
    best <- min(best, dt / k)
  }
  best
}
set.seed(8101)
mk <- function(n, ng) {
  d <- data.frame(x = rnorm(n), g = factor(sample(ng, n, TRUE)))
  d$y <- 2 * exp(0.3 * d$x) + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5)
  d
}
cases <- list(
  n2e3 = list(f = bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
              d = mk(2000, 50)),
  n2e4 = list(f = bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
              d = mk(20000, 200)),
  n2e4_plain = list(f = bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE),
                    d = mk(20000, 200)))
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.1))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- suppressWarnings(frm(cs$f, data = cs$d, start = st))
  fr <- fit$frame
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, st)
  t_guard <- block(function() ns$nl_check_body_read(fr, nll, tpl))
  t_dom <- block(function() ns$nl_domain_start(fr, tpl, integer(0), FALSE))
  t_fit <- block(function() suppressWarnings(frm(cs$f, data = cs$d, start = st)))
  cat(sprintf("%-11s guard %.4f s  domain %.4f s  fit %.4f s  (guard + domain) / fit = %.3f\n",
              nm, t_guard, t_dom, t_fit, (t_guard + t_dom) / t_fit))
}
