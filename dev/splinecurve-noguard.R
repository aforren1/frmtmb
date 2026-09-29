# Lane splinecurve: see the difference-guard test FAIL without the guard.
#
# Copies the lane's frmtmb.spline source to a scratch directory, deletes
# the `has_extra && isTRUE(allow_new_levels)` refusal from
# sp_curve_parts(), loads that copy with pkgload over the lane's private
# library, and runs test-new-levels.R against it.
#   Rscript dev/splinecurve-noguard.R > dev/splinecurve-noguard.log
lib <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline"
scratch <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/splinecurve-check/noguard"
unlink(scratch, recursive = TRUE)
dir.create(scratch, recursive = TRUE)
file.copy(src, scratch, recursive = TRUE)
pkg <- file.path(scratch, "frmtmb.spline")
f <- file.path(pkg, "R", "curve-cov.R")
x <- readLines(f)
start <- grep("if [(]has_extra && isTRUE[(]allow_new_levels[)][)] [{]", x)
stopifnot(length(start) == 1L)
# the refusal ends at the first line that is exactly "  }" after it
end <- start + which(x[(start + 1L):length(x)] == "  }")[1L]
cat("removing lines", start, "to", end, "of curve-cov.R:\n")
cat(x[start:end], sep = "\n")
writeLines(x[-(start:end)], f)
stopifnot(!any(grepl("isTRUE[(]allow_new_levels[)]", readLines(f))))

suppressMessages(library(testthat))
pkgload::load_all(pkg, quiet = TRUE)
res <- test_file(file.path(pkg, "tests", "testthat", "test-new-levels.R"),
                 package = "frmtmb.spline", load_package = "none",
                 env = testthat::test_env("frmtmb.spline"),
                 reporter = "silent")
for (t in res) {
  for (ex in t$results) {
    if (inherits(ex, c("expectation_failure", "expectation_error"))) {
      cat("FAILED [", t$test, "] ", conditionMessage(ex), "\n", sep = "")
    }
  }
}
r <- as.data.frame(res)
cat("RESULT noguard test-new-levels.R pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped), "\n",
    sep = "")

# What the unguarded call returns, in numbers, on the vignette's data
# (set.seed(4)) with (1 | subject) in place of the fs term.
set.seed(4)
n_sub <- 20
n_rep <- 12
n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(
  subject = factor(sub),
  trial = rep(seq_len(n_sub * n_rep), each = n_t),
  t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
fit <- frmtmb::frm(frmtmb::bf(v ~ s(t, k = 12) + (1 | subject)),
                   family = gaussian(), data = d)
lev <- c(levels(d$subject), "A", "B")
gA <- data.frame(t = seq(0, 1, length.out = 5),
                 subject = factor("A", levels = lev))
gB <- transform(gA, subject = factor("B", levels = lev))
cv <- frm_curve(fit, newdata = gA, contrast = gB, re_formula = NULL,
                allow_new_levels = TRUE, simultaneous = FALSE)
ev <- frmtmb::frm_lp_basis(fit, newdata = gA, re_formula = NULL,
                           allow_new_levels = TRUE)$extra_var
cat("unguarded difference of two unseen levels: .se =",
    format(cv$.se), "\n  the right se is at least sqrt(2 * extra_var) =",
    format(sqrt(2 * ev[1])), "\n")
