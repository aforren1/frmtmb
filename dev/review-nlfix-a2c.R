# Area 2: what the guard and the domain check cost on large nl models, in one
# process, block timed. Seed 6401. Also the unstr/fcor route, whose objective on
# doubles is dense.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
block <- function(f, min_t = 0.5, rounds = 3) {
  best <- Inf
  for (r in seq_len(rounds)) {
    k <- 1L
    repeat {
      t0 <- proc.time()[["elapsed"]]
      for (i in seq_len(k)) f()
      dt <- proc.time()[["elapsed"]] - t0
      if (dt > min_t) break
      k <- k * 2L
    }
    best <- min(best, dt / k)
  }
  best
}
one <- function(lab, form, data, fam = gaussian(), start = NULL, data2 = NULL) {
  t_fit <- system.time(fit <- suppressWarnings(suppressMessages(frm(form, family = fam, data = data, start = start, data2 = data2))))[["elapsed"]]
  fr <- fit$frame
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, start)
  t_g <- block(function() ns$nl_check_body_read(fr, nll, tpl))
  t_d <- block(function() ns$nl_domain_start(fr, tpl, integer(0), FALSE))
  t_b <- block(function() ns$build_objective(fr))
  cat(sprintf("%-26s n=%-7d fit %8.3f s  guard %7.4f s  domain %7.4f s (build %.4f)  (g+d)/fit = %.3f\n", lab, nrow(data), t_fit, t_g, t_d, t_b, (t_g + t_d) / t_fit))
}
set.seed(6401)
mk <- function(n, ng) {
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(ng, n, TRUE)))
  d$y <- 2 * exp(0.3 * d$x) + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5)
  d$t <- ave(seq_len(n), d$g, FUN = seq_along)
  d <- d[order(d$g, d$t), ]
  d
}
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.1))
one("plain 2e3", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), mk(2000, 50), start = st)
one("plain 2e5", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), mk(2e5, 200), start = st)
one("plain 1e6", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), mk(1e6, 200), start = st)
one("RE 2e4", bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE), mk(2e4, 200), start = st)
one("RE 2e5", bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE), mk(2e5, 500), start = st)
one("no start, linear body", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), mk(2e5, 200))
d <- mk(4000, 100)
one("ar cov 4e3", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), d, start = st)
one("ar cond 4e3", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), d, start = st)
d3 <- mk(300, 3)
M <- 0.3^abs(outer(seq_len(300), seq_len(300), "-"))
one("fcor dense 300", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE, autocor = ~ fcor(M)), d3, start = st, data2 = list(M = M))
d4 <- mk(3000, 30); d4$yc <- rpois(3000, exp(0.3 + 0.2 * d4$x))
one("poisson 3e3 RE", bf(yc ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), d4, fam = poisson())
d5 <- mk(30000, 30); d5$yo <- cut(0.8 * d5$x + rlogis(30000), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
one("cumulative 3e4", bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), d5, fam = cumulative())
