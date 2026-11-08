# Review of lane mfxse, second short check (S1, s1): hypotheses() naming
# group-level rows, a function hypothesis of the coefficients alone, and the
# s1 arming on seen levels. Data seed 7411 (mre_data()).
#   Rscript dev/review-mfxse-sc2.R
.libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib", "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
quiet <- function(e) suppressWarnings(suppressMessages(e))
try_ <- function(e) tryCatch(quiet(e), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 150)))
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
u0 <- rnorm(J, 0, 0.8); u1 <- rnorm(J, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x + 0.3 * (d$f == "b") + u0[d$g]))
fit <- quiet(frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d))
se_v <- sqrt(diag(vcov(fit)))
cat("sqrt(diag(vcov(fit))):", sprintf("%.6f", se_v), "\n")
show <- function(lab, r) {
  if (is.character(r)) cat(sprintf("%-45s %s\n", lab, r)) else
    cat(sprintf("%-45s rows %d  est %s  se %s\n", lab, nrow(r),
                paste(sprintf("%.5f", utils::head(r$estimate, 4)), collapse = " "),
                paste(sprintf("%.6f", utils::head(r$std.error, 4)), collapse = " ")))
}
cat("\n== hypotheses naming group-level rows\n")
show("'b4 = 0' (b4 past the 3 coefficients)", try_(hypotheses(fit, "b4 = 0")))
nm <- names(get_coef(fit))[4]
show(paste0("'`", nm, "` = 0'"), try_(hypotheses(fit, paste0("`", nm, "` = 0"))))
show("~ pairwise", try_(hypotheses(fit, hypothesis = ~ pairwise)))
show("~ reference", try_(hypotheses(fit, hypothesis = ~ reference)))
show("'b1 + b2 = 0'", try_(hypotheses(fit, "b1 + b2 = 0")))
cat("\n== a function hypothesis of the coefficients alone\n")
r1 <- try_(hypotheses(fit, hypothesis = function(x) coef(x)[1:2]))
show("function(x) coef(x)[1:2]", r1)
cat("  class of coef(fit):", class(coef(fit))[1], "\n")
r2 <- try_(hypotheses(fit, hypothesis = function(x) data.frame(term = c("i", "x"), estimate = fixef(x)[1:2, "Estimate"])))
show("function(x) fixef(x)[1:2, 'Estimate']", r2)
r3 <- try_(hypotheses(fit, hypothesis = function(x) data.frame(term = c("i", "x"), estimate = get_coef(x)[1:2])))
show("function(x) get_coef(x)[1:2]", r3)
if (!is.character(r2)) cat("  ratio to vcov():", sprintf("%.6f", r2$std.error / se_v[1:2]), "\n")
if (!is.character(r3)) cat("  ratio to vcov():", sprintf("%.6f", r3$std.error / se_v[1:2]), "\n")
cat("\n== s1: arming under vcov = ~ g\n")
wc <- function(e) {
  k <- 0; ws <- character(0)
  withCallingHandlers(e, warning = function(w) {
    if (inherits(w, "frmtmb_interop_se_part")) { k <<- k + 1; ws <<- c(ws, substr(conditionMessage(w), 1, 70)) }
    invokeRestart("muffleWarning")
  })
  paste(k, paste(ws, collapse = " | "))
}
f1 <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
nds <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "4", "5")))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
cat("predictions seen, ~ g              ", wc(predictions(f1, newdata = nds, vcov = ~ g)), "\n")
cat("avg_slopes seen, ~ g               ", wc(avg_slopes(f1, newdata = nds, variables = "x", vcov = ~ g)), "\n")
cat("avg_comparisons seen, ~ g, anl     ", wc(avg_comparisons(f1, newdata = nds, variables = "f", vcov = ~ g, allow_new_levels = TRUE)), "\n")
cat("avg_predictions all data, ~ g      ", wc(avg_predictions(f1, vcov = ~ g)), "\n")
cat("full vcov, seen                    ", wc(predictions(f1, newdata = nds, vcov = vcov(f1, full = TRUE))), "\n")
cat("~ g at new levels                  ", wc(predictions(f1, newdata = ndn, vcov = ~ g, allow_new_levels = TRUE)), "\n")
cat("then predictions seen, vcov = FALSE", wc(predictions(f1, newdata = nds, vcov = FALSE)), "\n")
cat("then avg_slopes seen, ~ g          ", wc(avg_slopes(f1, newdata = nds, variables = "x", vcov = ~ g)), "\n")
cat("hypotheses(f1, vcov = ~ g) then predictions new, vcov = FALSE ",
    { invisible(quiet(hypotheses(f1, vcov = ~ g))); wc(predictions(f1, newdata = ndn, vcov = FALSE, allow_new_levels = TRUE)) }, "\n")
cat("  armed after hypotheses(vcov = ~ g):", isTRUE(f1$cache[["mfx_rows_armed"]]), "\n")
