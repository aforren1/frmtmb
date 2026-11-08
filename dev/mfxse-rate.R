# The off-rate of marginaleffects' default standard errors on mixed fits,
# the review's construction (dev/review-predread-r2-mfxrate.R) on the lane's
# library setup: for each design, the delta-method SE against the SD of 150
# draws of the whole joint covariance (fixed effects, variance parameters and
# group effects), each read through the same quantity. Ratio = marginaleffects
# / reference; "material" is outside 0.8 to 1.25. fitted()'s Est.Error is
# kept beside the 12 levels as a second reference.
#   MFX_ARM=before|after Rscript dev/mfxse-rate.R <family> <out.tsv>
# Seeds 7501 to 7503 per design; draw seed = data seed + 1000.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
a <- commandArgs(TRUE)
famname <- a[1]
outf <- a[2]
fams <- list(gaussian = gaussian(), poisson = poisson(),
             bernoulli = bernoulli(), lognormal = lognormal(),
             negbinomial = negbinomial(), binomial = binomial())
fam <- fams[[famname]]
gen <- function(eta, n) {
  switch(famname,
         gaussian = eta + rnorm(n, 0, 0.8),
         poisson = rpois(n, exp(eta)),
         bernoulli = rbinom(n, 1, plogis(eta)),
         lognormal = exp(eta + rnorm(n, 0, 0.5)),
         negbinomial = rnbinom(n, mu = exp(eta), size = 3),
         binomial = rbinom(n, 10, plogis(eta)))
}
designs <- expand.grid(J = c(8L, 30L), m = c(10L, 40L), sdg = c(0.2, 0.8),
                       slope = c(FALSE, TRUE), stringsAsFactors = FALSE)
designs <- designs[!(designs$J == 8L & designs$m == 40L &
                       designs$sdg == 0.2), ]
rows <- list()
for (dd in seq_len(nrow(designs))) for (s in 7501:7503) {
  J <- designs$J[dd]
  m <- designs$m[dd]
  sdg <- designs$sdg[dd]
  sl <- designs$slope[dd]
  set.seed(s)
  n <- J * m
  d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(J)), each = m)),
                  x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
  gw <- rnorm(J)
  d$w <- gw[d$g]
  u0 <- rnorm(J, 0, sdg)[d$g]
  u1 <- if (sl) rnorm(J, 0, sdg / 2)[d$g] else 0
  eta <- 0.4 + (0.5 + u1) * d$x + 0.4 * (d$f == "b") + 0.3 * d$w + u0
  d$y <- gen(eta, n)
  d$tr <- 10L
  rhs <- if (sl) "x + f + w + (1 + x | g)" else "x + f + w + (1 | g)"
  fml <- stats::as.formula(paste(if (famname == "binomial") {
    "y | trials(tr)"
  } else "y", "~", rhs))
  fit <- tryCatch(quiet(frm(bf(fml), family = fam, data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat("FITERR", famname, J, m, sdg, sl, s, conditionMessage(fit), "\n")
    next
  }
  jc <- frmtmb:::get_joint_cov(fit)
  V <- as.matrix(jc$V)
  nm <- jc$names
  if (any(!is.finite(V))) {
    cat("NONFINITE", famname, J, m, sdg, sl, s, "\n")
    next
  }
  ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
  Rt <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), length(ee$values))
  idx <- split(seq_along(nm), factor(nm, levels = unique(nm)))
  center <- unlist(fit$estimates[names(idx)], use.names = FALSE)
  at <- function(p) {
    f <- fit
    for (k in names(idx)) f$estimates[[k]][] <- p[idx[[k]]]
    f$cache <- new.env(parent = emptyenv())
    f
  }
  rs <- function(f, nd) {
    as.numeric(frm_linpred(f, newdata = nd, type = "response"))
  }
  ndp <- d[seq(1, n, length.out = 12), ]
  dh <- d
  dl <- d
  h <- 1e-4
  dh$x <- d$x + h
  dl$x <- d$x - h
  da <- d
  da$f <- factor("a", levels = levels(d$f))
  db <- d
  db$f <- factor("b", levels = levels(d$f))
  dw1 <- d
  dw0 <- d
  dw1$w <- d$w + 1
  qty <- function(f) {
    c(pred = rs(f, ndp), avg_pred = mean(rs(f, d)),
      slope_x = mean((rs(f, dh) - rs(f, dl)) / (2 * h)),
      cmp_f = mean(rs(f, db) - rs(f, da)), cmp_w = mean(rs(f, dw1) - rs(f, dw0)))
  }
  set.seed(s + 1000L)
  M <- 150L
  draws <- vapply(seq_len(M), function(i) {
    qty(at(center + as.numeric(Rt %*% rnorm(length(center)))))
  }, numeric(16L))
  ref <- apply(draws, 1, sd)
  ws <- character(0)
  me <- function(e) {
    tryCatch(withCallingHandlers(e, warning = function(w) {
      ws <<- c(ws, class(w)[1])
      invokeRestart("muffleWarning")
    }, message = function(mm) invokeRestart("muffleMessage")),
    error = function(err) NULL)
  }
  t0 <- proc.time()[["elapsed"]]
  pp <- me(marginaleffects::predictions(fit, newdata = ndp))
  ap <- me(marginaleffects::avg_predictions(fit))
  sx <- me(marginaleffects::avg_slopes(fit, variables = "x"))
  cf <- me(marginaleffects::avg_comparisons(fit, variables = "f"))
  cw <- me(marginaleffects::avg_comparisons(fit, variables = list(w = 1)))
  secs <- proc.time()[["elapsed"]] - t0
  fse <- tryCatch(quiet(fitted(fit, newdata = ndp))[, "Est.Error"],
                  error = function(e) rep(NA_real_, 12L))
  out <- data.frame(family = famname, J = J, m = m, sdg = sdg, slope = sl,
                    seed = s,
                    qty = c(rep("pred", 12L), "avg_pred", "slope_x", "cmp_f",
                            "cmp_w"),
                    me_se = c(if (is.null(pp)) rep(NA, 12L) else pp$std.error,
                              if (is.null(ap)) NA else ap$std.error,
                              if (is.null(sx)) NA else sx$std.error,
                              if (is.null(cf)) NA else cf$std.error,
                              if (is.null(cw)) NA else cw$std.error),
                    ref = ref, est = rowMeans(draws),
                    fitted_se = c(unname(fse), rep(NA, 4L)),
                    warn = paste(unique(ws), collapse = ","),
                    conv = fit$opt$convergence, secs = secs)
  out$ratio <- out$me_se / out$ref
  rows[[length(rows) + 1L]] <- out
  cat(famname, J, m, sdg, sl, s, " pred ratio",
      signif(range(out$ratio[1:12]), 3), " avg", signif(out$ratio[13:16], 3),
      " secs", secs, " warn", out$warn[1], "\n")
}
res <- do.call(rbind, rows)
write.table(res, outf, sep = "\t", row.names = FALSE, quote = FALSE)
cat("DONE", famname, nrow(res) / 16, "fits\n")
