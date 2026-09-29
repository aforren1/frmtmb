# Reviewer, lane splinecurve, claim 3: the new contrast refusal.
#   Rscript splinecurve-rev-04-contrast.R > splinecurve-rev-log/04-contrast.txt
# Runs each case on the lane's installed build (guarded), then loads a
# copy of the lane source with the refusal deleted (pkgload, in a
# scratch directory outside the tree) and reruns them.
# Data: the vignette simulation, set.seed(4), plus a within-subject
# two-level factor `cond` drawn with set.seed(11).
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
set.seed(4)
n_sub <- 20; n_rep <- 12; n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
set.seed(11)
d$cond <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
d$v <- d$v + 0.05 * (d$cond == "b")

fit_re <- frm(bf(v ~ s(t, k = 12) + (1 | subject)), family = gaussian(),
              data = d)
fit_rc <- frm(bf(v ~ cond + s(t, k = 12) + (1 | subject)),
              family = gaussian(), data = d)
fit_fs <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
              family = gaussian(), data = d)
lev <- c(levels(d$subject), "A", "B")
g5 <- data.frame(t = seq(0, 1, length.out = 5))
gA <- transform(g5, subject = factor("A", levels = lev))
gB <- transform(g5, subject = factor("B", levels = lev))
gAa <- transform(gA, cond = factor("a", levels = levels(d$cond)))
gAb <- transform(gA, cond = factor("b", levels = levels(d$cond)))
gS3 <- transform(g5, subject = factor("3", levels = lev))

cases <- function(tag) {
  r <- function(expr) tryCatch(expr, error = function(e) e)
  fmt <- function(label, x) {
    if (inherits(x, "error")) {
      cat(tag, label, ": REFUSED:", substr(conditionMessage(x), 1, 80), "\n")
    } else {
      cat(tag, label, ": est", format(x$.estimate, digits = 5), "\n",
          strrep(" ", nchar(tag) + nchar(label) + 2), "se ",
          format(x$.se, digits = 5), "\n")
    }
    invisible(x)
  }
  fmt("1 two unseen (1|g) levels A-B, NULL, TRUE",
      r(frm_curve(fit_re, newdata = gA, contrast = gB, re_formula = NULL,
                  allow_new_levels = TRUE, simultaneous = FALSE)))
  fmt("2 same unseen level A, cond b-a, NULL, TRUE",
      r(frm_curve(fit_rc, newdata = gAb, contrast = gAa, re_formula = NULL,
                  allow_new_levels = TRUE, simultaneous = FALSE)))
  fmt("2' same, re_formula = NA (the refusal's advice)",
      r(frm_curve(fit_rc, newdata = gAb, contrast = gAa, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)))
  fmt("3 two unseen fs levels A-B, TRUE",
      r(frm_curve(fit_fs, newdata = gA, contrast = gB,
                  allow_new_levels = TRUE, simultaneous = FALSE)))
  fmt("3' same, simultaneous = TRUE",
      r(frm_curve(fit_fs, newdata = gA, contrast = gB,
                  allow_new_levels = TRUE, simultaneous = TRUE)))
  x <- fmt("4 seen fs level 3 minus unseen fs level A, TRUE",
           r(frm_curve(fit_fs, newdata = gS3, contrast = gA,
                       allow_new_levels = TRUE, simultaneous = FALSE)))
  if (!inherits(x, "error")) {
    la <- frm_lp_basis(fit_fs, newdata = gS3, re_formula = NA)
    lb <- frm_lp_basis(fit_fs, newdata = gA, re_formula = NA,
                       allow_new_levels = TRUE)
    Dm <- as.matrix(la$A) - as.matrix(lb$A)
    se_h <- sqrt(diag(Dm %*% la$V %*% t(Dm)))
    e3 <- as.numeric(frm_linpred(fit_fs, newdata = gS3, re_formula = NA))
    eA <- as.numeric(frm_linpred(fit_fs, newdata = gA, re_formula = NA,
                                 allow_new_levels = TRUE))
    cat(tag, "4 hand (A3 - Anew) V (A3 - Anew)': max rel diff in se",
        format(max(abs(x$.se / se_h - 1)), digits = 3),
        " est vs linpred(3) - linpred(new): max abs",
        format(max(abs(x$.estimate - (e3 - eA))), digits = 3), "\n")
  }
}
cat("sd(subject), (1|subject) fit:",
    format(VarCorr(fit_re)$subject$sd[1, 1], digits = 6),
    " sqrt(2) * sd:", format(sqrt(2) * VarCorr(fit_re)$subject$sd[1, 1],
                               digits = 6), "\n")
cat("sd(subject), cond fit:",
    format(VarCorr(fit_rc)$subject$sd[1, 1], digits = 6), "\n\n")
cases("GUARDED  ")

# the lane's source with the refusal deleted
src <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline"
scratch <- file.path(Sys.getenv("TEMP"), "splinecurve-rev-noguard")
unlink(scratch, recursive = TRUE)
dir.create(scratch, recursive = TRUE)
file.copy(src, scratch, recursive = TRUE)
f <- file.path(scratch, "frmtmb.spline", "R", "curve-cov.R")
x <- readLines(f)
start <- grep("if [(]has_extra && isTRUE[(]allow_new_levels[)][)] [{]", x)
stopifnot(length(start) == 1L)
end <- start + which(x[(start + 1L):length(x)] == "  }")[1L]
writeLines(x[-(start:end)], f)
stopifnot(!any(grepl("isTRUE[(]allow_new_levels[)]", readLines(f))))
cat("\nremoved lines", start, "to", end, "of curve-cov.R\n")
pkgload::load_all(file.path(scratch, "frmtmb.spline"), quiet = TRUE,
                  export_all = FALSE)
cat("noguard frm_curve from:", environmentName(environment(frm_curve)), "\n")
cases("UNGUARDED")
