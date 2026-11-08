# Review of lane mfxse, attacks 4 and 5: user vcov = and the row warnings.
#   Rscript dev/review-mfxse-vcovwarn.R
# Data seed 7411 (mre_data()), 8150 (the exact gp design), 9301 (a
# heteroskedastic gaussian (1 | g) for the sandwich).
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lane, base, user))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
# run e, return the value and the frmtmb_interop_se_part warnings it raised
wcount <- function(e) {
  ws <- character(0)
  v <- withCallingHandlers(e, warning = function(w) {
    if (inherits(w, "frmtmb_interop_se_part")) ws <<- c(ws, substr(conditionMessage(w), 1, 80))
    invokeRestart("muffleWarning")
  })
  list(v = v, n = length(ws), w = ws)
}
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
u0 <- rnorm(J, 0, 0.8); u1 <- rnorm(J, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x + 0.3 * (d$f == "b") + u0[d$g]))
fit <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))

cat("\n== attack 5: the row warnings\n")
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
nds <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "4", "5")))
show <- function(lab, r) cat(sprintf("%-55s warnings %d %s\n", lab, r$n,
                                     paste(unique(r$w), collapse = " | ")))
show("predictions, 2 new-level rows", wcount(predictions(fit, newdata = ndn, allow_new_levels = TRUE)))
show("avg_predictions, 2 new-level rows", wcount(avg_predictions(fit, newdata = ndn, allow_new_levels = TRUE)))
show("avg_slopes(x), 2 new-level rows", wcount(avg_slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE)))
show("avg_comparisons(f), 2 new-level rows", wcount(avg_comparisons(fit, newdata = ndn, variables = "f", allow_new_levels = TRUE)))
show("predictions, seen levels, allow_new_levels = TRUE", wcount(predictions(fit, newdata = nds, allow_new_levels = TRUE)))
show("predictions, new levels, re_formula = NA", wcount(predictions(fit, newdata = ndn, re_formula = NA, allow_new_levels = TRUE)))
show("predictions, new levels, vcov = FALSE", wcount(predictions(fit, newdata = ndn, vcov = FALSE, allow_new_levels = TRUE)))
show("predictions, new levels, sample_new_levels = old_levels",
     wcount(tryCatch(predictions(fit, newdata = ndn, allow_new_levels = TRUE,
                                 sample_new_levels = "old_levels"),
                     error = function(e) {cat("  error:", conditionMessage(e), "\n"); NULL})))
# the armed flag after a call that forms a covariance and predicts nothing
invisible(wcount(hypotheses(fit, joint = TRUE)))
cat("armed after hypotheses(joint = TRUE):", isTRUE(fit$cache[["mfx_rows_armed"]]), "\n")
show("then predictions(vcov = FALSE, new levels)", wcount(predictions(fit, newdata = ndn, vcov = FALSE, allow_new_levels = TRUE)))
invisible(get_vcov(fit))
cat("armed after get_vcov(fit):", isTRUE(fit$cache[["mfx_rows_armed"]]), "\n")
show("then predictions(vcov = FALSE, new levels)", wcount(predictions(fit, newdata = ndn, vcov = FALSE, allow_new_levels = TRUE)))
# two warnings in one call: a new level and an unseen gp position
set.seed(8150)
mg <- 90
dg <- data.frame(x = runif(mg, 0, 6), g = factor(sample(1:6, mg, TRUE)))
dg$y <- sin(dg$x) + rnorm(6, 0, 0.5)[dg$g] + rnorm(mg, sd = 0.3)
fg <- quiet(frm(y ~ gp(x) + (1 | g), data = dg))
ndg <- data.frame(x = c(dg$x[1], 3.33, 9), g = factor(c("1", "new", "2")))
show("gp + (1|g): one new level, one far position", wcount(predictions(fg, newdata = ndg, allow_new_levels = TRUE)))
ft <- quiet(fitted(fg, newdata = ndg, allow_new_levels = TRUE))
pg <- quiet(predictions(fg, newdata = ndg, allow_new_levels = TRUE))
cat("  mfx/fitted:", sprintf("%.4f", pg$std.error / ft[, "Est.Error"]), "\n")
# the kriging share at points between the data: how far from the data a
# row must be before the 1e-4 share is crossed, and what it moves
xs <- sort(dg$x)
gaps <- diff(xs)
mid <- xs[-1] - gaps / 2
ndm <- data.frame(x = c(mid[order(-gaps)][1:3], max(xs) + c(0.05, 0.2, 0.5)), g = factor("1", levels = levels(dg$g)))
fm <- quiet(fitted(fg, newdata = ndm))
pm <- wcount(predictions(fg, newdata = ndm))
lb <- frm_lp_basis(fg, newdata = ndm)
th <- fg$estimates$theta
bk <- Filter(function(b) b$covstruct == "gp", fg$frame$re_blocks)[[1]]
v0 <- exp(th[bk$theta_idx[1]])^2
cat("  rows: largest 3 gaps' midpoints (gap", sprintf("%.3f", sort(gaps, TRUE)[1:3]), ") and max(x) + 0.05, 0.2, 0.5\n")
cat("  kriging share:", sprintf("%.1e", lb$extra_var / v0), "\n")
cat("  mfx/fitted:  ", sprintf("%.4f", pm$v$std.error / fm[, "Est.Error"]), "\n")
cat("  warned rows:", paste(pm$w, collapse = " | "), "\n")

cat("\n== attack 4: vcov =\n")
set.seed(9301)
J <- 30L; m <- 10L; n <- J * m
dh <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
dh$y <- 1 + 0.5 * dh$x + rnorm(J, 0, 0.7)[dh$g] + rnorm(n, 0, exp(0.6 * dh$x) * 0.5)
dh$yok <- 1 + 0.5 * dh$x + rnorm(J, 0, 0.7)[dh$g] + rnorm(n, 0, 0.5)
for (resp in c("yok", "y")) {
  fh <- quiet(frm(stats::as.formula(paste(resp, "~ x + (1 | g)")), data = dh))
  nd <- dh[c(1, 101, 201), ]
  p0 <- predictions(fh, newdata = nd)
  ph <- tryCatch(quiet(predictions(fh, newdata = nd, vcov = "HC0")), error = function(e) { cat("  HC0 refused:", substr(conditionMessage(e), 1, 90), "
"); list(std.error = NA) })
  pc <- quiet(predictions(fh, newdata = nd, vcov = ~ g))
  a0 <- avg_predictions(fh); ah <- tryCatch(quiet(avg_predictions(fh, vcov = "HC0")), error = function(e) list(std.error = NA)); ac <- quiet(avg_predictions(fh, vcov = ~ g))
  cat(sprintf("%s: level SE default %s | HC0/default %s | CR1/default %s\n", resp,
              paste(sprintf("%.4f", p0$std.error), collapse = " "),
              paste(sprintf("%.3f", ph$std.error / p0$std.error), collapse = " "),
              paste(sprintf("%.3f", pc$std.error / p0$std.error), collapse = " ")))
  cat(sprintf("   avg SE default %.4f HC0 %.4f CR1 %.4f\n", a0$std.error, ah$std.error, ac$std.error))
  # split of a level's variance: inner Qrr^-1 (model-based, not robust)
  # against the carried part J Vf J'
  Q <- ns$joint_precision(fh); r <- fh$obj$env$random
  Wrr <- as.matrix(Matrix::solve(Q[r, r]))
  jc <- ns$get_joint_cov(fh)
  bb <- which(jc$names == "b")
  cat(sprintf("   b block: mean diag inner %.4f, mean diag total (model) %.4f\n",
              mean(diag(Wrr)), mean(diag(as.matrix(jc$V)[bb, bb]))))
}
# a matrix over the coefficients alone: the warning, and what a carry of
# that matrix through the coefficients' rows of J would give instead
V0 <- get_vcov(fit)[1:3, 1:3]
r1 <- wcount(predictions(fit, newdata = nds, vcov = V0))
cat(sprintf("\ncoef-only matrix: warnings %d; SE %s; default SE %s\n", r1$n,
            paste(sprintf("%.4f", r1$v$std.error), collapse = " "),
            paste(sprintf("%.4f", predictions(fit, newdata = nds)$std.error), collapse = " ")))
# REML and a sandwich
fr <- quiet(frm(yok ~ x + (1 | g), data = dh, REML = TRUE))
rr <- wcount(tryCatch(predictions(fr, newdata = dh[1:3, ], vcov = ~ g),
                      error = function(e) conditionMessage(e)))
cat("REML vcov = ~ g:", if (is.character(rr$v)) rr$v else sprintf("%.4f", rr$v$std.error), " warnings", rr$n, "\n")
