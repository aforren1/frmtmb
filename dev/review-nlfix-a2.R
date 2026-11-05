# Area 2: the guard. False refusals, false passes, ordering, cost. Seed 6201.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
res <- function(lab, expr) {
  r <- tt(expr)
  if (iserr(r)) cat(sprintf("%-34s ERROR %s\n", lab, substr(r$msg, 1, 260)))
  else cat(sprintf("%-34s ll=%s conv=%s beta=%s %s\n", lab, format(as.numeric(logLik(r)), digits = 9),
                   r$opt$convergence, paste(format(r$estimates$beta, digits = 5), collapse = ","),
                   if (length(attr(r, "w"))) paste0("WARN[", substr(attr(r, "w")[1], 1, 120), "]") else ""))
  invisible(r)
}
set.seed(6201)
n <- 200
d <- data.frame(x = rnorm(n), x2 = runif(n, 10, 20), g = factor(rep(1:10, each = 20)))
d$y <- 0.5 * d$x + rnorm(n)
d$yb <- rbinom(n, 1, 0.3 + 0.4 * (d$x2 > 15))
cat("--- bodies with a flat region around the start\n")
# a change point past the data, one coefficient: step(x2 - c) at c = 0, .37, -.53 is 1 on all rows
d$ys <- ifelse(d$x2 > 15, 1, 0) + rnorm(n, 0, 0.3)
res("step(x2 - c), one coef", frm(bf(ys ~ step(x2 - c), c ~ 1, nl = TRUE), data = d))
res("step(x2 - c), start 15", frm(bf(ys ~ step(x2 - c), c ~ 1, nl = TRUE), data = d, start = list(beta = c("c_(Intercept)" = 15))))
res("a * step(x2 - c)", frm(bf(ys ~ a * step(x2 - c), a + c ~ 1, nl = TRUE), data = d))
res("pmax(x2 - c, 0)", frm(bf(ys ~ pmax(x2 - c, 0), c ~ 1, nl = TRUE), data = d))
res("pmin(x2, c), flat at start", frm(bf(ys ~ pmin(x2 - 8, c), c ~ 1, nl = TRUE), data = d, start = list(beta = c("c_(Intercept)" = 100))))
res("ifelse in body", frm(bf(ys ~ ifelse(x2 > 15, a, 0), a ~ 1, nl = TRUE), data = d))
cat("--- saturated bodies: bitwise equal, body is read\n")
d$yt <- exp(-0.002 * d$x2 * 100) + rnorm(n, 0, 0.1)
d$xh <- d$x2 * 1000
res("exp(-k xh) start k=1", frm(bf(yt ~ exp(-k * xh), k ~ 1, nl = TRUE), data = d, start = list(beta = c("k_(Intercept)" = 1))))
res("a tiny multiplier a*1e-300", frm(bf(y ~ a * x * 1e-300, a ~ 1, nl = TRUE), data = d))
d$xt <- d$x * 1e-20
res("tiny x, linear", frm(y ~ 0 + xt, data = d))
res("tiny x, a * xt", frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d))
cat("--- a null of the chosen direction\n")
res("2 * a - b", frm(bf(y ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = d))
res("2 * a - b + 0 * x, prior", frm(bf(y ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = d,
                           prior = c(set_prior("normal(0, 1)", nlpar = "a"), set_prior("normal(0, 1)", nlpar = "b"))))
res("(2 * a - b) * x", frm(bf(y ~ (2 * a - b) * x, a + b ~ 1, nl = TRUE), data = d))
res("a - 2 * b + c", frm(bf(y ~ a - 2 * b + c, a + b + c ~ 1, nl = TRUE), data = d))
res("a*b - b*a (really absent)", frm(bf(y ~ a * b - b * a, a + b ~ 1, nl = TRUE), data = d))
cat("--- a coefficient at a zero multiplier, partly\n")
d$x0 <- 0
res("a * x0", frm(bf(y ~ a * x0, a ~ 1, nl = TRUE), data = d))
res("a * x0, linear", frm(y ~ 0 + x0, data = d))
res("b + a * x0", frm(bf(y ~ b + a * x0, a + b ~ 1, nl = TRUE), data = d))
cat("--- a body flat at its start: flat product at zero\n")
res("a * b * x", frm(bf(y ~ a * b * x, a + b ~ 1, nl = TRUE), data = d))
cat("--- fixed and mapped\n")
res("constant prior on a", frm(bf(y ~ a * x, a ~ 1, nl = TRUE), data = d, prior = set_prior("constant(0.5)", nlpar = "a")))
res("constant prior, all", frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d, prior = c(set_prior("constant(0)", nlpar = "a"), set_prior("constant(0.5)", nlpar = "b"))))
res("a ~ 0 + offset", frm(bf(y ~ a * x, a ~ 0 + offset(rep(0.5, 200)), nl = TRUE), data = d))
res("bound coefficient", frm(bf(y ~ a * x, a ~ 1, nl = TRUE), data = d, prior = set_prior("", nlpar = "a", lb = 0.2, ub = 0.9)))
cat("--- body read only under a branch of a custom family\n")
lp_branch <- frmtmb_family("rvbranch", dpars = c("mu", "sigma"), links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    # the mean reaches only the rows with y above 2, about 2 in 100 rows
    m <- dpars$mu * as.numeric(y > 2)
    RTMB::dnorm(y, m, dpars$sigma, log = TRUE)
  })
res("custom, mu at y > 2 only", frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = lp_branch, data = d))
lp_branch0 <- frmtmb_family("rvbranch0", dpars = c("mu", "sigma"), links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    m <- dpars$mu * 0
    RTMB::dnorm(y, m, dpars$sigma, log = TRUE)
  })
res("custom, mu * 0", frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = lp_branch0, data = d))
lp_int <- frmtmb_family("rvhalf", dpars = c("mu", "sigma"), links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    # reads a constant-across-rows part of mu only: mean(mu)
    mm <- sum(dpars$mu) / length(y)
    RTMB::dnorm(y, mm, dpars$sigma, log = TRUE)
  })
res("custom, reads mean(mu) only", frm(bf(y ~ b1 * (x - mean(x)), b1 ~ 1, nl = TRUE), family = lp_int, data = d))
cat("--- a dpar that the family never reads, by a second route\n")
res("nlf(sigma) body", frm(bf(y ~ x, nlf(sigma ~ s0 + s1 * x), s0 + s1 ~ 1), data = d))
res("hurdle all-zero y", {
  dz <- d; dz$yz <- 0L
  frm(bf(yz ~ exp(b0 + b1 * x), b0 + b1 ~ 1, nl = TRUE), family = hurdle_poisson(), data = dz)
})
res("hurdle all-zero y linear", {
  dz <- d; dz$yz <- 0L
  frm(yz ~ x, family = hurdle_poisson(), data = dz)
})
cat("--- sampler route and the guard\n")
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  suppressPackageStartupMessages(library(frmtmb.sample))
  r <- tt(frm_sample(bf(y ~ step(x2 - c), c ~ 1, nl = TRUE), data = d, chains = 1, iter = 20, refresh = 0,
                     prior = set_prior("normal(0, 10)", nlpar = "c"), sample_prior = "only"))
  show("frm_sample flat body, prior only", if (iserr(r)) r else 1)
  r <- tt(frm_sample(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = lp_branch0, data = d, chains = 1, iter = 20, refresh = 0))
  show("frm_sample body the family skips", if (iserr(r)) r else 1)
}
cat("--- ordering against a start out of domain\n")
res("log(a) default start", frm(bf(y ~ log(a) * x, a ~ 1, nl = TRUE), data = d))
res("custom skip + log(a)", frm(bf(y ~ log(a) * x, a ~ 1, nl = TRUE), family = lp_branch0, data = d))
res("custom skip, gamma-like", frm(bf(y ~ 1 / b0, b0 ~ 1, nl = TRUE), family = lp_branch0, data = d))
cat("--- par_template, start and control path\n")
res("control silent", frm(bf(y ~ a * x, a ~ 1, nl = TRUE), data = d, control = frmtmb_control(check_se = FALSE)))
