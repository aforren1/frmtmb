# Review of lane slopetrig, check 4: back-transforms on fits the ratio
# rule newly engages. Per case: D (default), T (autoscale = TRUE), Z (the
# same model on the centered, standardized column, x = m + sx * xz).
# D against T: identical() on every accessor. D against Z mapped back:
# the largest relative difference per accessor (sd and cor through the
# covariance map A = [[1, -m/sx], [0, 1/sx]]; the slope by 1/sx; the
# intercept by the slope times m/sx).
# Usage: Rscript slopetrig-review-back.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
mk <- function(seed, s, c0, fam = "gau") {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- stats::rnorm(300)
  eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    stats::rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  y <- if (fam == "gau") eta + stats::rnorm(300) else stats::rpois(300, exp(eta))
  data.frame(g = g, x = c0 + s * z, y = y)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
rel <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(NA)
  max(abs(a - b)) / max(abs(b))
}
cases <- list(
  list("year s3 seed2 int", 2, 3, 2000, y ~ 1 + (1 + x | g), "gau"),
  list("year s3 seed2 fix", 2, 3, 2000, y ~ x + (1 + x | g), "gau"),
  list("5+0.1z seed1 int", 1, 0.1, 5, y ~ 1 + (1 + x | g), "gau"),
  list("5+0.1z seed3 fix", 3, 0.1, 5, y ~ x + (1 + x | g), "gau"),
  list("age 40+10z fix", 1, 10, 40, y ~ x + (1 + x | g), "gau"),
  list("c0=10 ratio2.2 fix", 4, 10 / sqrt(2.2^2 - 1), 10, y ~ x + (1 + x | g), "gau"),
  list("Days2-11-like 1+3z fix", 5, 3 / 1, 6.5, y ~ x + (1 + x | g), "gau"),
  list("pois 5+0.1z fix", 2, 0.1, 5, y ~ x + (1 + x | g), "pois"))
for (cs in cases) {
  d <- mk(cs[[2]], cs[[3]], cs[[4]], cs[[6]])
  fam <- if (cs[[6]] == "gau") gaussian() else poisson()
  fo <- cs[[5]]
  m <- mean(d$x); sx <- stats::sd(d$x)
  dz <- d; dz$x <- (d$x - m) / sx
  fd <- q(frm(fo, data = d, family = fam))
  ft <- q(frm(fo, data = d, family = fam,
              control = frmtmb_control(autoscale = TRUE)))
  fz <- q(frm(fo, data = dz, family = fam))
  cat("\n==", cs[[1]], ": engaged", isTRUE(fd$autoscaled),
      "; ratio", signif(sqrt(mean(d$x^2)) / sx, 4),
      "; ll D - Z", signif(as.numeric(logLik(fd)) - as.numeric(logLik(fz)), 3),
      "\n")
  nd <- data.frame(g = d$g[c(1, 16, 31)], x = m + sx * c(-1.5, 0, 2))
  ndz <- nd; ndz$x <- (nd$x - m) / sx
  acc <- list(
    par = function(f, d) f$opt$par,
    fixef = function(f, d) fixef(f),
    vcov = function(f, d) vcov(f),
    VarCorr = function(f, d) q(VarCorr(f)),
    ranef = function(f, d) ranef(f),
    coef = function(f, d) coef(f),
    fitted = function(f, d) fitted(f),
    pred_new = function(f, d) { set.seed(9); q(predict(f, newdata = d)) },
    fitted_new_nore = function(f, d) q(fitted(f, newdata = d, re_formula = NA)),
    draws = function(f, d) { set.seed(9); q(predict(f, newdata = d, re_formula = NA, probs = c(0.05, 0.95))) },
    hypothesis = function(f, d) { set.seed(9); q(hypothesis(f, "Intercept > 0")) },
    ce = function(f, d) { set.seed(9); q(conditional_effects(f, "x")) })
  for (nm in names(acc)) {
    a <- tryCatch(acc[[nm]](fd, nd), error = function(e) e)
    b <- tryCatch(acc[[nm]](ft, nd), error = function(e) e)
    same <- if (inherits(a, "error") || inherits(b, "error")) {
      paste("ERROR", conditionMessage(if (inherits(a, "error")) a else b))
    } else {
      # compare(): drop the call/environment parts that differ by object
      identical(unclass(a), unclass(b)) ||
        isTRUE(all.equal(a, b, tolerance = 0, check.attributes = FALSE))
    }
    cat(sprintf("  D vs T %-16s %s\n", nm, format(same)))
  }
  # D against Z mapped back
  vz <- q(VarCorr(fz))$g; vd <- q(VarCorr(fd))$g
  sdz <- vz$sd[, "Estimate"]; corz <- vz$cor[1, "Estimate", 2]
  Sz <- diag(sdz) %*% matrix(c(1, corz, corz, 1), 2) %*% diag(sdz)
  A <- matrix(c(1, 0, -m / sx, 1 / sx), 2)
  Sr <- A %*% Sz %*% t(A)
  sdr <- sqrt(diag(Sr)); corr <- Sr[1, 2] / prod(sdr)
  cat(sprintf("  D vs Z sd %.3g ; cor %.3g (cor D %.10f, mapped %.10f) ; 1+cor rel %.3g\n",
              rel(vd$sd[, "Estimate"], sdr), rel(vd$cor[1, "Estimate", 2], corr),
              vd$cor[1, "Estimate", 2], corr,
              rel(1 + abs(vd$cor[1, "Estimate", 2]), 1 + abs(corr))))
  # VarCorr Est.Error of the sds: delta method on Z is not mapped here;
  # report D's Est.Error and whether it is finite
  cat("  D VarCorr Est.Error sd:", signif(vd$sd[, "Est.Error"], 4),
      "cor:", signif(vd$cor[1, "Est.Error", 2], 4), "\n")
  rd <- as.matrix(ranef(fd)$g); rz <- as.matrix(ranef(fz)$g)
  cat(sprintf("  D vs Z ranef slope %.3g ; ranef intercept %.3g\n",
              rel(rd[, "x"] * sx, rz[, "x"]),
              rel(rd[, "Intercept"], rz[, "Intercept"] - rz[, "x"] * m / sx)))
  fe_d <- fixef(fd); fe_z <- fixef(fz)
  if ("x" %in% rownames(fe_d)) {
    Vz <- vcov(fz)[c("Intercept", "x"), c("Intercept", "x")]
    B <- matrix(c(1, 0, -m / sx, 1 / sx), 2)
    br <- B %*% fe_z[c("Intercept", "x"), "Estimate"]
    Vr <- B %*% Vz %*% t(B)
    cat(sprintf("  D vs Z fixef est %.3g ; fixef SE %.3g\n",
                rel(fe_d[c("Intercept", "x"), "Estimate"], br),
                rel(fe_d[c("Intercept", "x"), "Est.Error"], sqrt(diag(Vr)))))
  } else {
    cat(sprintf("  D vs Z fixef est %.3g ; fixef SE %.3g\n",
                rel(fe_d[, "Estimate"], fe_z[, "Estimate"]),
                rel(fe_d[, "Est.Error"], fe_z[, "Est.Error"])))
  }
  cat(sprintf("  D vs Z fitted est %.3g ; fitted SE %.3g ; newdata (no re) est %.3g SE %.3g\n",
              rel(fitted(fd)[, "Estimate"], fitted(fz)[, "Estimate"]),
              rel(fitted(fd)[, "Est.Error"], fitted(fz)[, "Est.Error"]),
              rel(q(fitted(fd, newdata = nd, re_formula = NA))[, "Estimate"],
                  q(fitted(fz, newdata = ndz, re_formula = NA))[, "Estimate"]),
              rel(q(fitted(fd, newdata = nd, re_formula = NA))[, "Est.Error"],
                  q(fitted(fz, newdata = ndz, re_formula = NA))[, "Est.Error"])))
  ced <- tryCatch(q(conditional_effects(fd, "x"))[[1]], error = function(e) NULL)
  cez <- tryCatch(q(conditional_effects(fz, "x"))[[1]], error = function(e) NULL)
  if (!is.null(ced) && !is.null(cez) && nrow(ced) == nrow(cez)) {
    cat(sprintf("  D vs Z conditional_effects x grid %.3g ; estimate %.3g ; lower %.3g\n",
                rel(ced$x, m + sx * cez$x), rel(ced$estimate__, cez$estimate__),
                rel(ced$lower__, cez$lower__)))
  }
  hd <- tryCatch(q(hypothesis(fd, "Intercept > 0"))$hypothesis, error = function(e) NULL)
  hz <- tryCatch(q(hypothesis(fz, "Intercept > 0"))$hypothesis, error = function(e) NULL)
  if (!is.null(hd)) cat("  hypothesis D:", format(unlist(hd[1, 2:4]), digits = 6),
                        "; Z (Intercept at the mean of x unless fixed x absent):",
                        format(unlist(hz[1, 2:4]), digits = 6), "\n")
  # draws of the population-level fitted values, 4000 each: column means
  # and sds, D against Z (Monte Carlo noise about 2 percent on the sd)
  set.seed(11)
  dd <- tryCatch(q(fitted(fd, newdata = nd, re_formula = NA, summary = FALSE,
                          ndraws = 4000)), error = function(e) NULL)
  set.seed(12)
  dz2 <- tryCatch(q(fitted(fz, newdata = ndz, re_formula = NA, summary = FALSE,
                           ndraws = 4000)), error = function(e) NULL)
  if (!is.null(dd) && !is.null(dz2)) {
    cat(sprintf("  D vs Z fitted draws: mean %.3g ; sd %.3g\n",
                rel(colMeans(dd), colMeans(dz2)),
                rel(apply(dd, 2, sd), apply(dz2, 2, sd))))
  }
}
cat("DONE\n")
