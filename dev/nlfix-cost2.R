# Minor 1, re-measured on the review's cases (dev/review-nlfix-a2c.R, seed 6401),
# with the guard as frm() now runs it: the taped gradient first, the objective
# only for what is left, nl_domain_start() on the objective already built, and
# the saddle's look at the body.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
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
  t_fit <- system.time(fit <- suppressWarnings(suppressMessages(
    frm(form, family = fam, data = data, start = start, data2 = data2))))[["elapsed"]]
  fr <- fit$frame
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, start)
  # the taped objective at its start, as the guard reads it
  ob <- suppressWarnings(frm(form, family = fam, data = data, start = start,
                             data2 = data2, dry_run = "objective"))$obj
  t_g <- block(function() {
    g <- ns$nl_guard_gradient(ob, fr, NULL)
    ns$nl_check_body_read(fr, nll, tpl, g)
  })
  t_d <- block(function() ns$nl_domain_start(fr, tpl, integer(0), FALSE, nll))
  t_s <- block(function() ns$nl_saddle_find(fit))
  cat(sprintf(paste0("%-24s n=%-7d fit %8.3f s  guard %7.4f s  domain %7.4f s  ",
                     "saddle look %7.4f s  (all)/fit = %.3f\n"), lab, nrow(data), t_fit,
              t_g, t_d, t_s, (t_g + t_d + t_s) / t_fit))
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
one("RE 2e4", bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE), mk(2e4, 200), start = st)
one("no start, linear body", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), mk(2e5, 200))
d <- mk(4000, 100)
one("ar cov 4e3", bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE,
                     autocor = ~ ar(time = t, gr = g, cov = TRUE)), d, start = st)
d5 <- mk(30000, 30)
d5$yo <- cut(0.8 * d5$x + rlogis(30000), c(-Inf, -1, 0, 1, Inf), labels = 1:4,
             ordered_result = TRUE)
one("cumulative 3e4", bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), d5, fam = cumulative())
one("cumulative 3e4, saddle", bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), d5,
    fam = cumulative())
cat("\nthe restart's cost, where it runs: a saddle fit against a started fit, data of\n")
cat("dev/nlfix-saddle.R, seed 9202\n")
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf),
                   labels = 1:4))
f <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
t1 <- block(function() suppressWarnings(suppressMessages(frm(f, family = cumulative(), data = d))))
t2 <- block(function() suppressWarnings(frm(f, family = cumulative(), data = d,
  start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.1)))))
cat(sprintf("saddle fit %.3f s, started fit %.3f s, ratio %.2f\n", t1, t2, t1 / t2))
