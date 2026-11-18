# Short check of lane squant fix round 1, items 4 and 5: fits behind
# frm_compat("quantile_elf") rows marked "works" (data and settings not
# the lane's), and the doc's standard-error routes.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-compat.R > dev/reviews/squant-short-compat.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ct <- frm_compat("quantile_elf")
cat("frm_compat('quantile_elf') (core only):", nrow(ct), "rows;",
    paste(names(table(ct$status)), table(ct$status), collapse = ", "), "\n")
other <- ifelse(ct$feature_a == "quantile_elf", ct$feature_b, ct$feature_a)
pick <- c("nl", "mixture", "us", "diag", "cs", "s()", "t2()", "gp_pred()",
          "mi()")
for (k in pick) {
  i <- which(other == k)
  cat(sprintf("  row %-10s %-12s %s\n", k,
              if (length(i)) ct$status[i[1]] else "(none)",
              if (length(i)) substr(ct$note[i[1]], 1, 150) else ""))
}
cat("\n")
run <- function(lab, expr) {
  ws <- character()
  ms <- character()
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 150))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, substr(conditionMessage(m), 1, 150))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error")) {
    cat(sprintf("%s: ERROR %s\n", lab, substr(conditionMessage(r), 1, 300)))
    return(invisible(NULL))
  }
  se <- tryCatch(sqrt(diag(vcov(r))), error = function(e) NA)
  cat(sprintf("%s: conv %s | logLik %.4f | warnings %d | messages %d | max SE %.3g | class %s | %.1f s\n",
              lab, r$opt$convergence, as.numeric(logLik(r)), length(ws),
              length(ms), max(se), paste(class(r), collapse = "/"), el))
  if (length(ws)) cat("   W:", paste(unique(ws), collapse = "\n   W: "), "\n")
  if (length(ms)) cat("   M:", paste(unique(ms), collapse = "\n   M: "), "\n")
  fe <- fixef(r)
  print(round(fe[seq_len(min(nrow(fe), 6)), 1:2, drop = FALSE], 4))
  invisible(r)
}
set.seed(20261117)
n <- 400
d <- data.frame(x = runif(n), x2 = runif(n),
                g = factor(sample(sprintf("g%02d", 1:25), n, TRUE)))
re0 <- rnorm(25, 0, 0.6)
re1 <- rnorm(25, 0, 0.5)
gi <- as.integer(d$g)
d$y <- 1 + 2 * d$x + re0[gi] + re1[gi] * d$x + (0.5 + d$x) * rnorm(n)
d$ys <- sin(2 * pi * d$x) + 0.5 * d$x2^2 + 0.4 * rnorm(n)
d$xp <- 0.2 + 3 * d$x
d$yn <- 2 * d$xp^0.7 + 0.3 * rnorm(n)
z <- rbinom(n, 1, 0.5)
d$ym <- ifelse(z == 1, -3 + rnorm(n), 3 + rnorm(n, 0, 0.7))

run("nl (y = a xp^b, q 0.75)",
    frm(bf(yn ~ a * xp^b, a ~ 1, b ~ 1, nl = TRUE, quantile = 0.75),
        family = quantile_elf(), data = d,
        prior = c(prior(normal(2, 2), nlpar = "a"),
                  prior(normal(1, 1), nlpar = "b"))))
run("mixture (quantiles 0.25, 0.75)",
    frm(bf(ym ~ 1, mu1 ~ 1, mu2 ~ 1, quantile1 = 0.25, quantile2 = 0.75),
        family = mixture(quantile_elf(), quantile_elf()), data = d))
run("us: (1 + x | g), q 0.5, smoothing 1",
    frm(bf(y ~ x + (1 + x | g), quantile = 0.5),
        family = quantile_elf(smoothing = 1), data = d))
run("diag: diag(1 + x | g), q 0.75, smoothing 1",
    frm(bf(y ~ x + diag(1 + x | g), quantile = 0.75),
        family = quantile_elf(smoothing = 1), data = d))
run("cs: cs(1 + x | g), q 0.5, smoothing 1",
    frm(bf(y ~ x + cs(1 + x | g), quantile = 0.5),
        family = quantile_elf(smoothing = 1), data = d))
run("s(x) + x2, q 0.9",
    frm(bf(ys ~ s(x) + x2, quantile = 0.9), family = quantile_elf(),
        data = d))
run("t2(x, x2), q 0.25",
    frm(bf(ys ~ t2(x, x2), quantile = 0.25), family = quantile_elf(),
        data = d))
run("gp(x) + x2, q 0.5",
    frm(bf(ys ~ gp(x) + x2, quantile = 0.5), family = quantile_elf(),
        data = d[1:150, ]))

cat("\n== docs: the standard-error routes ==\n")
set.seed(1)
nn <- 300
dd <- data.frame(x = runif(nn), row = seq_len(nn))
dd$y <- 1 + 2 * dd$x + (0.5 + dd$x) * rnorm(nn)
fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(), data = dd)
V <- vcov_cluster(fit, ~ row, "CR0", full = TRUE)
s0 <- summary(fit)
s1 <- summary(fit, vcov = V)
print(s0$fixed[, 1:2])
print(s1$fixed[, 1:2])
print(confint(fit, vcov = V))
vc <- vcov(fit, cluster = ~ row)
Vd <- vcov_cluster(fit, ~ row)
cat("vcov(fit, cluster = ~ row) SEs:", sprintf("%.6f", sqrt(diag(vc))), "\n")
cat("vcov_cluster(fit, ~ row) default type, its fixed block SEs:",
    sprintf("%.6f", sqrt(diag(Vd))[seq_len(nrow(vc))]), "\n")
cat("CR0 full SEs (coef):", sprintf("%.6f", sqrt(diag(V))[seq_len(nrow(vc))]),
    "\n")
set.seed(3)
dg <- data.frame(x = runif(120), g = factor(rep(1:30, each = 4)))
dg$y <- 1 + dg$x + rnorm(30, 0, 0.5)[dg$g] + rnorm(120)
fq <- run("quadrature (1 | g), smoothing 1",
          frm(bf(y ~ x + (1 | g), quantile = 0.5),
              family = quantile_elf(smoothing = 1), data = dg,
              quadrature = TRUE))
if (!is.null(fq)) {
  r <- tryCatch(vcov_cluster(fq, ~ g, "CR0", full = TRUE),
                error = function(e) e)
  cat("vcov_cluster on the quadrature fit:",
      if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
        "returned a matrix", "\n")
}
fl <- run("Laplace (1 | g), smoothing 1",
          frm(bf(y ~ x + (1 | g), quantile = 0.5),
              family = quantile_elf(smoothing = 1), data = dg))
if (!is.null(fl)) {
  r <- tryCatch(vcov_cluster(fl, ~ g, "CR0", full = TRUE),
                error = function(e) e)
  cat("vcov_cluster(fit, ~ g, 'CR0', full = TRUE) on the Laplace fit:",
      if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
        sprintf("SEs %s", paste(sprintf("%.4f", sqrt(diag(r))[1:2]),
                                collapse = " ")), "\n")
}
cat("DONE\n")
