# Code-8 stops that are NOT at the maximum, to show the gain test still
# warns there. Two kinds:
# - forced: nlminb's own false-convergence test made coarse
#   (optCtrl xf.tol), with restarts = 0 so frmtmb keeps the stop; models
#   from GLMMs, location-scale and nonlinear fits;
# - natural: asym_laplace and zero_inflated_asym_laplace (the check-loss
#   kink; test-v14.R's data), whose stops are short of the exact
#   quantile-regression maximum.
# Each is judged against the same model's default fit (dll < 0 = short).
#   CODE8_ARM=before|after Rscript dev/code8-short.R <label> [seeds] [xf.tols]
a <- commandArgs(TRUE)
lab <- a[1]
# optional: seeds and xf.tol values of the forced part, comma lists
seeds <- if (length(a) > 1) as.integer(strsplit(a[2], ",")[[1]]) else 1:4
xfs <- if (length(a) > 2) as.numeric(strsplit(a[3], ",")[[1]]) else
  c(1e-6, 1e-4, 1e-2)
source("dev/code8-common.R")
rows <- list()
rec <- function(case, r, ref_ll) {
  si <- stop_info(r)
  rows[[length(rows) + 1L]] <<- cbind(
    data.frame(case = case), si,
    dll = if (r$err) NA else ll_of(r$f) - ref_ll,
    w_code = sum(grepl("frmtmb_conv_code", r$wc)),
    m_atmax = sum(grepl("frmtmb_converged_at_max", r$mc)),
    w_other = sum(!grepl("frmtmb_conv_code", r$wc)),
    conds = cond_txt(r))
}
models <- list(
  lmm = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:12, each = 15)), x = rnorm(180))
    u <- rnorm(12, 0, 0.8); v <- rnorm(12, 0, 0.4)
    d$y <- 1 + 0.5 * d$x + u[d$g] + v[d$g] * d$x + rnorm(180)
    list(f = bf(y ~ x + (1 + x | g)), fam = gaussian(), d = d)
  },
  pglmm = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200))
    d$y <- rpois(200, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.7)[d$g]))
    list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d)
  },
  lsstud = function(s) {
    set.seed(s)
    d <- data.frame(x = rnorm(300))
    d$y <- 1 + 0.5 * d$x + exp(0.2 + 0.4 * d$x) * rt(300, 5)
    list(f = bf(y ~ x, sigma ~ x), fam = student(), d = d)
  },
  nlexp = function(s) {
    set.seed(s)
    d <- data.frame(z = runif(200, 0, 5))
    d$y <- 3 * exp(-0.7 * d$z) + 0.5 + rnorm(200, 0, 0.2)
    list(f = bf(y ~ A * exp(-k * z) + c0, A + k + c0 ~ 1, nl = TRUE),
         fam = gaussian(), d = d,
         start = list(beta = c(1, 0.2, 0)))
  },
  bglmm = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:30, each = 8)), x = rnorm(240))
    d$y <- rbinom(240, 1, plogis(-0.2 + 0.8 * d$x + rnorm(30, 0, 1)[d$g]))
    list(f = bf(y ~ x + (1 | g)), fam = bernoulli(), d = d)
  })
fitm <- function(m, ctl = NULL) {
  args <- list(m$f, family = m$fam, data = m$d)
  if (!is.null(m$start)) args$start <- m$start
  if (!is.null(ctl)) args$control <- ctl
  do.call(rfit, args)
}
for (nm in names(models)) for (s in seeds) {
  m <- models[[nm]](s)
  ref <- fitm(m)
  ref_ll <- if (ref$err) NA else ll_of(ref$f)
  for (xf in xfs) {
    r <- fitm(m, frmtmb_control(restarts = 0,
                                optCtrl = list(iter.max = 1000,
                                               eval.max = 1000,
                                               xf.tol = xf)))
    rec(sprintf("%s s%d xf.tol=%g", nm, s, xf), r, ref_ll)
  }
}
# natural: the check-loss kink, against quantreg's exact maximum
ald_ll <- function(d, b, tau, zi = NULL) {
  # profile sigma out of the asymmetric Laplace at fixed b (closed form)
  r <- d$y - b[1] - b[2] * d$x
  rho <- r * (tau - (r < 0))
  s <- mean(rho)
  sum(log(tau * (1 - tau)) - log(s) - rho / s)
}
for (s in c(13, 14, 15, 16)) {
  set.seed(s); n <- 500; x <- rnorm(n)
  d <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n) * (1 + 0.4 * x^2))
  r <- rfit(bf(y ~ x, quantile = 0.25), family = asym_laplace(), data = d)
  ref <- if (requireNamespace("quantreg", quietly = TRUE)) {
    ald_ll(d, coef(quantreg::rq(y ~ x, tau = 0.25, data = d)), 0.25)
  } else NA
  rec(sprintf("asym_laplace s%d", s), r, ref)
}
X <- do.call(rbind, rows)
dir.create("dev/code8-out", showWarnings = FALSE)
write.table(X, sprintf("dev/code8-out/short-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(table(code = X$code, useNA = "ifany"))
c8 <- X[!is.na(X$code) & grepl("false", X$msg), ]
cat("fits", nrow(X), "false conv", nrow(c8), "\n")
print(c8[, c("case", "gain", "se_fin", "n_lost", "dll", "w_code", "m_atmax",
             "w_other")], row.names = FALSE)
cat("CONDS\n")
for (i in seq_len(nrow(c8))) cat(c8$case[i], ":", substr(c8$conds[i], 1, 300),
                                 "\n")
cat("DONE\n")
