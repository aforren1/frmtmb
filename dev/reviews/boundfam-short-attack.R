# Short check of lane boundfam's fix round: other sources of par_units
# (B1) and the separation rule for probability dpars (M1).
#   Rscript dev/reviews/boundfam-short-attack.R <lane|base> [b1|m1]
arm <- commandArgs(trailingOnly = TRUE)[1]
part <- commandArgs(trailingOnly = TRUE)[-1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/reviews/boundfam-review-fp.R")
eval(parse(text = src[grep("^conds <- ", src):(grep("^ok <- ", src) - 1L)]))
F <- function(...) conds(frm(...))
ok <- function(p) !length(part) || p %in% part
pu <- function(r) cat("  par_units:", if (is.null(r$value$par_units)) "NULL" else
  format(signif(r$value$par_units, 3)), "\n")

if (ok("b1")) {
  set.seed(11)
  n <- 300
  d <- data.frame(x = rnorm(n), z = runif(n),
                  o = factor(sample(1:4, n, TRUE), ordered = TRUE))
  d$xs <- d$x * 1e-6
  # mo(): a fit with a maximum, and a negbinomial at its poisson limit
  d$y <- rpois(n, exp(0.5 + 0.2 * as.integer(d$o) + 0.3 * d$x))
  r <- F(y ~ mo(o) + x, family = poisson(), data = d); show("mo poisson", r); pu(r)
  r <- F(y ~ mo(o) + x, family = negbinomial(), data = d); show("mo negbinomial on poisson", r); pu(r)
  r <- F(y ~ mo(o) + xs, family = poisson(), data = d,
         control = frmtmb_control(autoscale = FALSE)); show("mo + xs autoscale FALSE", r); pu(r)
  # me()
  d$sdx <- 0.3
  d$xm <- d$x + rnorm(n, 0, 0.3)
  d$yg <- 1 + 0.5 * d$x + rnorm(n, 0, 0.5)
  r <- F(yg ~ me(xm, sdx), data = d); show("me gaussian", r); pu(r)
  # gp()
  dg <- d[1:80, ]
  dg$yg <- sin(3 * dg$z) + rnorm(80, 0, 0.3)
  r <- F(yg ~ gp(z), data = dg); show("gp gaussian", r); pu(r)
  dg$yp <- rpois(80, exp(0.5 + sin(3 * dg$z)))
  r <- F(yp ~ xs + gp(z), family = poisson(), data = dg,
         control = frmtmb_control(autoscale = FALSE)); show("gp + xs autoscale FALSE", r); pu(r)
  # cs() in an adjacent-category ordinal
  d$yo <- factor(cut(0.6 * d$x + rlogis(n), c(-Inf, -1, 0.5, Inf), labels = FALSE), ordered = TRUE)
  r <- F(yo ~ cs(x), family = acat(), data = d); show("acat cs(x)", r); pu(r)
  # nonlinear
  d$yn <- 2 * exp(0.4 * d$x) + rnorm(n, 0, 0.3)
  r <- F(bf(yn ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
         start = list(beta = c(1.5, 0.3))); show("nonlinear", r); pu(r)
  # autoscale on (default plans xs) with a smooth: at the maximum, and a
  # true limit beside it
  d$yp <- rpois(n, exp(0.5 + 0.3 * d$x + sin(2 * pi * d$z)))
  r <- F(yp ~ xs + s(z), family = poisson(), data = d); show("autoscale default xs + s(z)", r); pu(r)
  r <- F(yp ~ xs + s(z), family = negbinomial(), data = d); show("autoscale default xs + s(z) negbinomial", r); pu(r)
  r <- F(yp ~ x + s(z), family = negbinomial(), data = d,
         control = frmtmb_control(autoscale = FALSE)); show("x + s(z) negbinomial autoscale FALSE", r); pu(r)
  # a numeric column with 3 values, the top one on 3 rows (small spread)
  d$t3 <- c(rep(0, n - 20), rep(1, 17), rep(2, 3))
  d$yb <- rbinom(n, 1, plogis(-1 + 0.5 * d$t3))
  r <- F(yb ~ t3 + x, family = bernoulli(), data = d); show("t3 rare level", r); pu(r)
  ref("glm", glm(yb ~ t3 + x, binomial, d))
}

if (ok("m1")) {
  set.seed(21)
  n <- 300
  d <- data.frame(f = factor(rep(c("a", "b"), each = n / 2)), x = rnorm(n))
  d$xc <- d$x - mean(d$x)
  # group a all zero, group b plain poisson
  d$y1 <- ifelse(d$f == "a", 0L, rpois(n, 4))
  show("zi ~ 0 + f, a all zero", F(bf(y1 ~ 1, zi ~ 0 + f), family = zero_inflated_poisson(), data = d))
  show("zi ~ f, a all zero", F(bf(y1 ~ 1, zi ~ f), family = zero_inflated_poisson(), data = d))
  d$z1 <- as.integer(d$y1 == 0)
  show("bernoulli analogue z ~ 0 + f", F(z1 ~ 0 + f, family = bernoulli(), data = d))
  # no excess zeros anywhere, centered continuous zi predictor
  d$y2 <- rpois(n, 5)
  show("zi ~ xc, no excess zeros", F(bf(y2 ~ 1, zi ~ xc), family = zero_inflated_poisson(), data = d))
  # no zeros at all, hurdle with a centered predictor: hu -> 0 for every row
  d$y3 <- 1L + rpois(n, 3)
  show("hu ~ xc, no zeros", F(bf(y3 ~ 1, hu ~ xc), family = hurdle_poisson(), data = d))
  # mixture weight separated by the sign of x
  d$y4 <- ifelse(d$x < 0, rnorm(n, -3, 0.5), rnorm(n, 3, 0.5))
  show("mixture theta1 ~ x separating", F(bf(y4 ~ 1, theta1 ~ x), family = mixture(gaussian, gaussian),
                                          data = d, start = list(beta = c(-3, 3))))
  # zero_one_inflated_beta: zoi separated by x; coi separated by x
  m <- plogis(0.2 + 0.3 * d$x)
  yb <- rbeta(n, m * 10, (1 - m) * 10)
  d$y5 <- ifelse(d$x > 0 & runif(n) < 0.3, as.numeric(runif(n) < 0.5), yb)
  show("zoi ~ x separating", F(bf(y5 ~ x, zoi ~ x), family = zero_one_inflated_beta(), data = d))
  d$y6 <- ifelse(runif(n) < 0.3, as.numeric(d$x > 0), yb)
  show("coi ~ x separating", F(bf(y6 ~ x, coi ~ x), family = zero_one_inflated_beta(), data = d))
}
