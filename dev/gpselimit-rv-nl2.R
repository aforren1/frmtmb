# Lane gpselimit, fix round 1: a copy of dev/reviews/gpselimit-review-nl2.R
# whose output goes to dev/gpselimit-out/rv/ (the review's logs stay).
# Review of lane gpselimit, check 5 (variant 2: body a * (1 - exp(-b * z)),
# whose default start b = 0 is a saddle). A field in a nonlinear parameter. The
# white message is now a message, which nl_fit_off_saddle() does not hold
# and nl_fit_probed() says for a first fit that a restart or a sign probe
# replaced. Searches fits of y ~ a * exp(b * z), a ~ gp(x, k), b ~ 1 for
# a frmtmb_gp_basis condition beside a restart, probe or pole message, or
# more than one frmtmb_gp_basis condition, and prints whether the returned
# fit's block is white.
# Usage: Rscript gpselimit-review-nl.R <lane|base> [seeds]
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-rv-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
seeds <- if (length(a) >= 2) eval(parse(text = a[2])) else 1:4
rv_libs(arm)
ns <- asNamespace("frmtmb")
for (b0 in c(1.5, -0.8)) for (k in c(10, 20)) for (seed in seeds) {
  set.seed(52300 + seed)
  n <- 120L
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  z <- abs(z) * 2
  y <- (2 + fs) * (1 - exp(-b0 * z)) + 0.1 * stats::rnorm(n)
  d <- data.frame(x, z, y)
  t0 <- proc.time()[["elapsed"]]
  q <- rv_fit(frm(bf(y ~ a * (1 - exp(-b * z)), a ~ gp(x, k = k), b ~ 1,
                     nl = TRUE), data = d, family = gaussian()))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(q$fit, "error")) {
    cat("NL", arm, b0, k, seed, "ERROR", conditionMessage(q$fit), "\n")
    next
  }
  fit <- q$fit
  bk <- rv_hsgp(fit)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  w <- ns$gp_basis_white(bk, th)
  ngp <- sum(vapply(q$conds, function(c) "frmtmb_gp_basis" %in% c$class, NA))
  other <- any(vapply(q$conds, function(c) {
    any(grepl("nl_|restart|pole", c$class))
  }, NA))
  cat(sprintf(paste0("NL %s b0 %4.1f k %2d seed %d | ll %.6f conv %d | ",
                     "ls/sp %.3g white %d | n_gp_basis %d nl_cond %d | ",
                     "%.1fs | conds [%s]\n"),
              arm, b0, k, seed, as.numeric(logLik(fit)),
              fit$opt$convergence, exp(th[2]) / bk$gp_spacing, w, ngp,
              other, el, rv_cstr(q$conds)))
  if (ngp > 1 || (ngp > 0 && other) || (ngp > 0 && !w)) {
    for (c in q$conds) {
      cat("   ", c$type, paste(c$class[1:2], collapse = "/"), ":",
          substr(gsub("\n", " ", c$text), 1, 160), "\n")
    }
  }
}
