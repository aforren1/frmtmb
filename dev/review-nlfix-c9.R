# Fix round 1: the rest of the first review's findings, rerun. Seeds as named.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
ll <- function(r) if (iserr(r)) paste("ERR", substr(r$msg, 1, 170)) else paste0(format(as.numeric(logLik(r)), digits = 9), " conv=", r$opt$convergence, " beta=", paste(format(r$estimates$beta, digits = 4), collapse = ","), " w=", length(attr(r, "w")), " m=", length(attr(r, "m")))
set.seed(6201); n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)), t = rep(1:20, 10), z = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
d$xt <- d$x * 1e-20
cat("tiny x  linear", ll(tt(frm(y ~ 0 + xt, data = d))), "\n  body  ", ll(tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d))), "\n")
d$xs <- d$x * 1e-10; d$ybig <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
cat("large y linear", ll(tt(frm(ybig ~ 0 + xs, data = d))), "\n  body  ", ll(tt(frm(bf(ybig ~ b * xs, b ~ 1, nl = TRUE), data = d))), "\n")
d$xm <- d$x * 1e-5; d$yr <- 0.5 * d$x + rnorm(n, 0, 0.2)
cat("x 1e-5  linear", ll(tt(frm(yr ~ 0 + xm, data = d))), "\n  body  ", ll(tt(frm(bf(yr ~ a * xm, a ~ 1, nl = TRUE), data = d))), "\n")
cat("rellib arm is", arm, "\n")
if (arm == "before") quit("no")
F <- ns$nl_stan_funs
cat("hypot(1e200,1e200)", F$hypot(1e200, 1e200), " hypot(1e-200,1e-200)", F$hypot(1e-200, 1e-200), " inv_square(0) grad", RTMB::MakeTape(function(p) F$inv_square(p), 0)$jacobian(0), " inv_cloglog grad at 5", RTMB::MakeTape(function(p) F$inv_cloglog(p), 5)$jacobian(5), "\n")
d$inv_logit <- d$x; d$yg <- plogis(3 * d$x) + rnorm(n, 0, 0.2)
cat("column and call:", ll(tt(frm(bf(yg ~ inv_logit(a * inv_logit), a ~ 1, nl = TRUE), data = d, start = list(beta = c("a_(Intercept)" = 1))))), "\n"); d$inv_logit <- NULL
d$yy <- 1 + 0.5 * d$x + ave(rnorm(n, 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
cat("nl = TRUE ar  :", ll(tt(frm(bf(yy ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), data = d))), "\n")
cat("nlf(mu) + ac  :", ll(tt(frm(bf(yy ~ 1) + nlf(mu ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1) + acformula(~ ar(time = t, gr = g, cov = TRUE)), data = d))), "\n")
cat("bf(autocor=) + nlf:", ll(tt(frm(bf(yy ~ 1, autocor = ~ ar(time = t, gr = g, cov = TRUE)) + nlf(mu ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1), data = d))), "\n")
for (f in c("rdm", "lba", "bandit2arm_delta", "hmm", "cox", "royston_parmar", "wiener", "wiener_gng", "gddm")) {
  suppressPackageStartupMessages({ library(frmtmb.eam); library(frmtmb.learn); library(frmtmb.survival); library(frmtmb.latent) })
  cat(sprintf("%-18s %s\n", f, as.data.frame(frm_compat("nl", f))$status))
}
cat("habit.Rmd:", grep("nl = TRUE", readLines("vignettes/habit.Rmd"), value = TRUE)[1], "\n")
cat("mixture(Gamma, Gamma) default start:", ll(tt(frm(bf(y1 ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = mixture(Gamma(), Gamma()), data = transform(d, y1 = exp(0.3 * x + rnorm(n, 0, 0.4)))))), "\n")
