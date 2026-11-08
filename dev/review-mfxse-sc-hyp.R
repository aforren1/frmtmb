# Review of lane mfxse, short check of the fix round: M1 (interop_in_hypotheses()),
# M2 (the row check armed at the estimates), m7's premise (which dots reach
# get_coef()). Data seed 7411 (mre_data()).
#   Rscript dev/review-mfxse-sc-hyp.R
.libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib", "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
u0 <- rnorm(J, 0, 0.8); u1 <- rnorm(J, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x + 0.3 * (d$f == "b") + u0[d$g]))
fit <- quiet(frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d))
nd <- d[c(1, 50, 120), ]
jt <- function(h) sprintf("F %.4f df1 %d", h$statistic, as.integer(h$df1))
cat("\n== M1: the joint test by every route (rellib-r18: F 15.4386 df1 3)\n")
w1 <- function(mm) marginaleffects::hypotheses(mm, joint = TRUE)
w2 <- function(object) hypotheses(object, joint = TRUE)
hh <- marginaleffects::hypotheses
cat("attached              ", jt(quiet(hypotheses(fit, joint = TRUE))), "\n")
cat("namespace-qualified   ", jt(quiet(marginaleffects::hypotheses(fit, joint = TRUE))), "\n")
cat("model = by name       ", jt(quiet(hypotheses(model = fit, joint = TRUE))), "\n")
cat("wrapper, other name   ", jt(quiet(w1(fit))), "\n")
cat("wrapper 2             ", jt(quiet(w2(fit))), "\n")
cat("alias hh()            ", jt(quiet(hh(fit, joint = TRUE))), "\n")
cat("do.call               ", jt(quiet(do.call(marginaleffects::hypotheses, list(fit, joint = TRUE)))), "\n")
cat("joint first, model 2nd", jt(quiet(hypotheses(joint = TRUE, fit))), "\n")
cat("\n== M1: the group-level effects must stay in these (SE equal to the plain call)\n")
cmp <- function(lab, a, b) cat(sprintf("%-50s %s | plain %s\n", lab,
  paste(sprintf("%.5f", a), collapse = " "), paste(sprintf("%.5f", b), collapse = " ")))
ap <- avg_predictions(fit); as <- avg_slopes(fit, variables = "x")
pr <- predictions(fit, newdata = nd)
cmp("hypotheses(avg_predictions(fit))", quiet(hypotheses(avg_predictions(fit)))$std.error, ap$std.error)
cmp("hypotheses(avg_slopes(fit))", quiet(hypotheses(avg_slopes(fit, variables = "x")))$std.error, as$std.error)
a1 <- avg_slopes(fit, variables = "x")
cmp("hypotheses(<slopes object>)", quiet(hypotheses(a1))$std.error, as$std.error)
p12 <- predictions(fit, newdata = nd, hypothesis = "b1 = b2")
cmp("hypotheses(predictions(nd), 'b1 = b2')", quiet(hypotheses(predictions(fit, newdata = nd), "b1 = b2"))$std.error, p12$std.error)
hf <- quiet(hypotheses(fit, hypothesis = function(x) {
  data.frame(term = "avg", estimate = avg_predictions(x)$estimate)
}))
cmp("hypotheses(fit, function(x) avg_predictions(x))", hf$std.error, ap$std.error)
hf2 <- quiet(hypotheses(fit, hypothesis = function(x) {
  data.frame(term = c("r1", "r2", "r3"), estimate = predictions(x, newdata = nd)$estimate)
}))
cmp("hypotheses(fit, function(x) predictions(x, nd))", hf2$std.error, pr$std.error)
both <- function(mm) {
  h <- hypotheses(mm, joint = TRUE)
  a <- avg_predictions(mm)
  list(h = h, a = a)
}
r <- quiet(both(fit))
cmp("function: hypotheses() then avg_predictions()", r$a$std.error, ap$std.error)
cat("   its joint test:", jt(r$h), "\n")
cat("\n== M1: cost of the frame walk\n")
deep <- function(k) if (k == 0) {
  t0 <- proc.time()[["elapsed"]]; for (i in 1:2000) ns$interop_in_hypotheses()
  (proc.time()[["elapsed"]] - t0) / 2000
} else deep(k - 1)
for (k in c(0, 50, 200)) cat(sprintf("  depth +%d: %.1f us per call\n", k, 1e6 * deep(k)))
cnt <- 0
tr <- function(...) cnt <<- cnt + 1
trace("interop_in_hypotheses", tr, where = ns, print = FALSE)
invisible(quiet(avg_slopes(fit, variables = "x"))); n_sl <- cnt; cnt <- 0
invisible(quiet(predictions(fit, newdata = nd))); n_pr <- cnt; cnt <- 0
untrace("interop_in_hypotheses", where = ns)
t_sl <- min(replicate(3, system.time(quiet(avg_slopes(fit, variables = "x")))[["elapsed"]]))
cat(sprintf("  calls per avg_slopes() %d, per predictions() %d; avg_slopes() %.2f s\n", n_sl, n_pr, t_sl))

cat("\n== M2: frmtmb_interop_se_part warnings reaching the caller\n")
wc <- function(e) {
  k <- 0
  withCallingHandlers(e, warning = function(w) {
    if (inherits(w, "frmtmb_interop_se_part")) k <<- k + 1
    invokeRestart("muffleWarning")
  })
  k
}
f1 <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
nds <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "4", "5")))
cat("predictions new     ", wc(predictions(f1, newdata = ndn, allow_new_levels = TRUE)), "\n")
cat("avg_predictions new ", wc(avg_predictions(f1, newdata = ndn, allow_new_levels = TRUE)), "\n")
cat("avg_slopes new      ", wc(avg_slopes(f1, newdata = ndn, variables = "x", allow_new_levels = TRUE)), "\n")
cat("avg_comparisons new ", wc(avg_comparisons(f1, newdata = ndn, variables = "f", allow_new_levels = TRUE)), "\n")
cat("slopes new          ", wc(slopes(f1, newdata = ndn, variables = "x", allow_new_levels = TRUE)), "\n")
cat("comparisons new     ", wc(comparisons(f1, newdata = ndn, variables = "f", allow_new_levels = TRUE)), "\n")
cat("avg_slopes seen     ", wc(avg_slopes(f1, newdata = nds, variables = "x", allow_new_levels = TRUE)), "\n")
invisible(get_vcov(f1))
cat("get_vcov(), then avg_slopes seen:            ", wc(avg_slopes(f1, newdata = nds, variables = "x")), "\n")
invisible(get_vcov(f1))
cat("get_vcov(), then predictions new vcov=FALSE: ", wc(predictions(f1, newdata = ndn, vcov = FALSE, allow_new_levels = TRUE)), "\n")
cat("   still armed:", isTRUE(f1$cache[["mfx_rows_armed"]]), "\n")
cat("   then hypotheses(f1, 'b2 = 0'):           ", wc(hypotheses(f1, "b2 = 0")), "\n")
cat("   then predictions new, re_formula = NA:    ", wc(predictions(f1, newdata = ndn, re_formula = NA, allow_new_levels = TRUE)), "\n")
invisible(get_vcov(f1))
V0 <- get_vcov(f1)[1:3, 1:3]
cat("get_vcov(), then predictions new, vcov = coef matrix: ",
    wc(predictions(f1, newdata = ndn, vcov = V0, allow_new_levels = TRUE)), "(the held warning and the stale-armed row count)\n")
invisible(get_vcov(f1))
cat("get_vcov(), then inferences(simulation) new: ",
    wc(tryCatch(inferences(predictions(f1, newdata = ndn, allow_new_levels = TRUE, vcov = FALSE), method = "simulation", R = 50),
                error = function(e) {cat("[error:", substr(conditionMessage(e), 1, 60), "] "); NULL})), "\n")

cat("\n== m7's premise: the dots each method receives, by call (allow_new_levels = TRUE)\n")
seen <- list()
gc0 <- ns$get_coef.frmtmb_fit; gv0 <- ns$get_vcov.frmtmb_fit; sc0 <- ns$set_coef.frmtmb_fit
registerS3method("get_coef", "frmtmb_fit", function(model, ...) {
  seen$get_coef <<- unique(c(seen$get_coef, paste(sort(names(list(...))), collapse = ","))); gc0(model, ...)
}, envir = asNamespace("marginaleffects"))
registerS3method("get_vcov", "frmtmb_fit", function(model, vcov = NULL, ...) {
  seen$get_vcov <<- unique(c(seen$get_vcov, paste(sort(names(list(...))), collapse = ","))); gv0(model, vcov, ...)
}, envir = asNamespace("marginaleffects"))
registerS3method("set_coef", "frmtmb_fit", function(model, coefs, ...) {
  seen$set_coef <<- unique(c(seen$set_coef, paste(sort(names(list(...))), collapse = ","))); sc0(model, coefs, ...)
}, envir = asNamespace("marginaleffects"))
for (call in c("predictions", "avg_predictions", "avg_slopes", "avg_comparisons")) {
  seen <- list()
  invisible(quiet(switch(call,
    predictions = predictions(f1, newdata = ndn, allow_new_levels = TRUE),
    avg_predictions = avg_predictions(f1, newdata = ndn, allow_new_levels = TRUE),
    avg_slopes = avg_slopes(f1, newdata = ndn, variables = "x", allow_new_levels = TRUE),
    avg_comparisons = avg_comparisons(f1, newdata = ndn, variables = "f", allow_new_levels = TRUE))))
  cat(sprintf("%-16s get_coef [%s] get_vcov [%s] set_coef [%s]\n", call,
              paste(seen$get_coef, collapse = " | "), paste(seen$get_vcov, collapse = " | "),
              paste(seen$set_coef, collapse = " | ")))
}
