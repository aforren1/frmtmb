# Reviewer: gradients of the censored and truncated objective against
# Richardson finite differences (numDeriv) of the same objective, at the
# start and with every distributional intercept moved to link-scale
# extremes; then a fit and its SEs against a numeric Hessian.
# GROUP = 1, 2 or 3 picks the families.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
grp <- as.integer(Sys.getenv("GROUP", "1"))
set.seed(20261020 + grp)
n <- 150
x <- rnorm(n)
codes <- sample(c(-1, 0, 0, 0, 1, 2), n, TRUE)
mk <- function(y, lo, hi, step, discrete, trials = NULL) {
  d <- data.frame(x = x, y = y, cc = codes)
  top <- if (!is.null(trials)) trials else hi
  d$y2 <- if (discrete) pmin(d$y + step, top) else pmin(d$y + step, hi)
  # an interval needs y2 > y (discrete: y2 >= y); drop the code where not
  bad <- d$cc == 2 & (d$y2 <= d$y)
  d$cc[bad] <- 0
  if (!is.null(trials)) d$tr <- trials
  d
}
fams <- list(
  student = list(fam = student(), d = mk(1 + 0.5 * x + rt(n, 4), -Inf, Inf, 0.7, FALSE)),
  skew_normal = list(fam = skew_normal(), d = mk(1 + 0.5 * x + abs(rnorm(n)), -Inf, Inf, 0.7, FALSE)),
  exgaussian = list(fam = exgaussian(), d = mk(1 + 0.5 * x + rnorm(n, 0, 0.3) + rexp(n), -Inf, Inf, 0.7, FALSE)),
  shifted_lognormal = list(fam = shifted_lognormal(), d = mk(0.3 + rlnorm(n, 0.2 * x, 0.5), 0, Inf, 0.4, FALSE)),
  asym_laplace = list(fam = asym_laplace(), d = mk(1 + 0.5 * x + rnorm(n), -Inf, Inf, 0.7, FALSE)),
  von_mises = list(fam = von_mises(), d = mk(pmin(pmax(2 * atan(0.3 * x) + rnorm(n, 0, 0.6), -3), 2.5), -pi, pi, 0.4, FALSE)),
  Gamma = list(fam = Gamma(link = "log"), d = mk(rgamma(n, 3, 3 / exp(0.5 + 0.3 * x)), 0, Inf, 0.5, FALSE)),
  inverse.gaussian = list(fam = inverse.gaussian(link = "log"), d = mk(statmod::rinvgauss(n, exp(0.3 * x), 3), 0, Inf, 0.5, FALSE)),
  Beta = list(fam = Beta(), d = mk(rbeta(n, 2 * exp(0.3 * x), 3), 0, 0.999, 0.1, FALSE)),
  binomial = list(fam = binomial(), d = mk(rbinom(n, 12, plogis(0.3 * x)), 0, 12, 2, TRUE, 12)),
  beta_binomial = list(fam = beta_binomial(), d = mk(rbinom(n, 12, rbeta(n, 2, 3)), 0, 12, 2, TRUE, 12)),
  negbinomial = list(fam = negbinomial(), d = mk(rnbinom(n, mu = exp(1 + 0.3 * x), size = 2), 0, Inf, 2, TRUE)),
  geometric = list(fam = geometric(), d = mk(rnbinom(n, mu = exp(1 + 0.3 * x), size = 1), 0, Inf, 2, TRUE)),
  poisson = list(fam = poisson(), d = mk(rpois(n, exp(1 + 0.3 * x)), 0, Inf, 2, TRUE)),
  compois = list(fam = brmsfamily("com_poisson"), d = mk(rpois(n, exp(1 + 0.3 * x)), 0, Inf, 2, TRUE)),
  zero_inflated_poisson = list(fam = zero_inflated_poisson(), d = mk(rpois(n, exp(1 + 0.3 * x)) * rbinom(n, 1, 0.8), 0, Inf, 2, TRUE)),
  zero_inflated_negbinomial = list(fam = zero_inflated_negbinomial(), d = mk(rnbinom(n, mu = exp(1 + 0.3 * x), size = 2) * rbinom(n, 1, 0.8), 0, Inf, 2, TRUE)),
  zero_inflated_binomial = list(fam = zero_inflated_binomial(), d = mk(rbinom(n, 12, 0.4) * rbinom(n, 1, 0.8), 0, 12, 2, TRUE, 12)),
  zero_inflated_beta_binomial = list(fam = zero_inflated_beta_binomial(), d = mk(rbinom(n, 12, rbeta(n, 2, 3)) * rbinom(n, 1, 0.8), 0, 12, 2, TRUE, 12)),
  hurdle_poisson = list(fam = hurdle_poisson(), d = mk(ifelse(runif(n) < 0.2, 0, 1 + rpois(n, exp(1 + 0.3 * x))), 0, Inf, 2, TRUE)),
  hurdle_negbinomial = list(fam = hurdle_negbinomial(), d = mk(ifelse(runif(n) < 0.2, 0, 1 + rnbinom(n, mu = exp(1 + 0.3 * x), size = 2)), 0, Inf, 2, TRUE)),
  hurdle_gamma = list(fam = hurdle_gamma(), d = mk(ifelse(runif(n) < 0.2, 0, rgamma(n, 3, 3 / exp(0.5 + 0.3 * x))), 0, Inf, 0.5, FALSE)),
  hurdle_lognormal = list(fam = hurdle_lognormal(), d = mk(ifelse(runif(n) < 0.2, 0, rlnorm(n, 0.2 * x, 0.6)), 0, Inf, 0.5, FALSE)),
  zero_inflated_beta = list(fam = zero_inflated_beta(), d = mk(ifelse(runif(n) < 0.2, 0, rbeta(n, 2, 3)), 0, 0.999, 0.1, FALSE))
)
pick <- names(fams)[(seq_along(fams) - 1) %% 3 + 1 == grp]
gcmp <- function(ob, p) {
  f0 <- ob$fn(p)
  g <- as.numeric(ob$gr(p))
  fd <- tryCatch(numDeriv::grad(ob$fn, p, method = "Richardson"),
                 error = function(e) rep(NA_real_, length(p)))
  c(fn = f0, err = max(abs(g - fd) / pmax(1, abs(fd))), nfg = sum(!is.finite(g)),
    nffd = sum(!is.finite(fd)))
}
for (nm in pick) {
  fm <- fams[[nm]]
  d <- fm$d
  tri <- nm %in% c("binomial", "beta_binomial", "zero_inflated_binomial",
                   "zero_inflated_beta_binomial")
  dpn <- setdiff(frmtmb:::as_frmtmb_family(fm$fam)$dpars, "mu")
  rhs <- paste(c("y |", if (tri) "trials(tr) +", "cens(cc, y2) ~ x"), collapse = " ")
  fo <- do.call(bf, c(list(as.formula(rhs)),
                      lapply(dpn, function(dp) as.formula(paste(dp, "~ 1")))))
  cat(sprintf("== %s (codes: %s)\n", nm, paste(names(table(d$cc)), table(d$cc), sep = ":", collapse = " ")))
  ob <- tryCatch(frm(fo, family = fm$fam, data = d, dry_run = "objective")$obj,
                 error = function(e) { cat("  BUILD ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(ob)) next
  p0 <- ob$par
  r <- gcmp(ob, p0)
  cat(sprintf("  start            fn %-12.6g grad vs FD %.2e  nonfinite grad %d fd %d\n",
              r["fn"], r["err"], r["nfg"], r["nffd"]))
  jd <- which(names(p0) == "betad")
  for (j in seq_along(jd)) for (v in c(-10, -4, 4, 10)) {
    p <- p0
    p[jd[j]] <- v
    r <- gcmp(ob, p)
    flag <- if (is.finite(r["fn"]) && (r["nfg"] > 0 || r["err"] > 1e-4)) "  <-- CHECK" else ""
    cat(sprintf("  %-6s eta = %3g  fn %-12.6g grad vs FD %.2e  nonfinite grad %d fd %d%s\n",
                dpn[j], v, r["fn"], r["err"], r["nfg"], r["nffd"], flag))
  }
  # a fit and its SEs against a numeric Hessian of the objective
  fit <- tryCatch(withCallingHandlers(frm(fo, family = fm$fam, data = d),
                                      warning = function(w) {
                                        cat("  WARN:", substr(conditionMessage(w), 1, 80), "\n")
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) { cat("  FIT ERROR:", substr(conditionMessage(e), 1, 120), "\n"); NULL })
  if (is.null(fit)) next
  H <- numDeriv::hessian(fit$obj$fn, fit$opt$par)
  se_num <- sqrt(diag(solve(H)))
  se_fit <- sqrt(diag(as.matrix(vcov(fit))))
  cat(sprintf("  fit conv %d logLik %.6f | max |grad| %.2e\n", fit$opt$convergence,
              as.numeric(logLik(fit)), max(abs(fit$obj$gr(fit$opt$par)))))
  cat(sprintf("  SE fit: %s\n  SE num: %s\n  max |ratio - 1| %.2e\n",
              paste(signif(se_fit, 6), collapse = " "),
              paste(signif(se_num, 6), collapse = " "),
              if (length(se_fit) == length(se_num)) max(abs(se_fit / se_num - 1)) else NA))
}
