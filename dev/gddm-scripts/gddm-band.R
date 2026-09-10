# lane gddm, punch round 1: the review's BLOCKER 1 construction, rerun
# against the per-pair tolerance.
#
# The first shipped tolerance was 1e-8 times the COLUMN MAXIMUM, applied
# to every pair of entries in the column. On a column whose dynamic
# range exceeds 1e8 the band is therefore wider than the small entries
# themselves, and the guard accepts a model whose small rows differ by
# whole units. The review built one: an intertemporal-choice delay of
# 1, 2 and 3 seconds beside ten years.
#
# The tolerance is now
#     1e-8 * max(|a|, |b|)  +  1e-11 * max|finite column|
# so the first term is scale-free and the second is only the floor a
# near-zero poly() entry needs.
#
# Seed 909, the review's own. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

set.seed(909)
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
n <- 120L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$blk <- rep(1:4, each = 30L)
# 1, 2, 3 s in the short blocks; ten years in the long ones
d$delay <- c(rep(c(1, 2, 3), each = 10L), rep(c(1, 2, 3), each = 10L),
             rep(3.15e8, 30L), rep(3.15e8, 30L))
d$cond <- d$blk
cat("delay column: max", format(max(d$delay), digits = 3),
    " old tolerance (1e-8 * colmax)", format(1e-8 * max(d$delay)), "\n")
cat("within-condition spread per block:",
    paste(tapply(d$delay, d$cond, function(z) diff(range(z))),
          collapse = " "), "\n")

form <- bf(rt | vint(upper, cond) ~ delay, bs ~ 1, ndt ~ 1, bias = 0.5)
r <- tryCatch({
  frm(form, family = gddm(control = ctl), data = d, dry_run = "frame")
  "ACCEPTED, no error and no warning"
}, error = function(e) paste("refused:", conditionMessage(e)))
cat("\n1. does the guard fire?  ", substr(r, 1, 200), "\n")

# Centering the same column: under the old tolerance this flipped the
# answer, which is the sharpest statement of the defect.
d2 <- d
d2$delay <- d$delay - mean(d$delay)
r2 <- tryCatch({
  frm(bf(rt | vint(upper, cond) ~ delay, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = gddm(control = ctl), data = d2, dry_run = "frame")
  "ACCEPTED"
}, error = function(e) "refused")
cat("2. the same column centered:", r2,
    "  (the two must agree)\n")

# If the guard still accepts, the defect is still there: two data sets
# differing on 116 of 120 rows with a bitwise equal objective.
if (!grepl("^refused", r)) {
  fit <- frm(form, family = gddm(control = ctl), data = d)
  p <- fit$opt$par
  first <- vapply(split(seq_len(n), d$cond), min, integer(1))
  ob <- function(dat) frm(form, family = gddm(control = ctl),
                          data = dat, dry_run = "objective")$obj
  fa <- as.numeric(ob(d)$fn(p))
  dB <- d
  j <- setdiff(seq_len(n), first)
  dB$delay[j] <- dB$delay[j] + 1
  fb <- as.numeric(ob(dB)$fn(p))
  cat(sprintf("3. fn(A) %.17g\n   fn(B) %.17g   identical: %s\n",
              fa, fb, identical(fa, fb)))
  cat("   coefficients:",
      paste(names(unlist(fixef(fit))),
            sprintf("%.6g", unlist(fixef(fit))), sep = "=",
            collapse = "  "), "\n")
}

# The tolerance the two smallest entries of that column now get.
cat("\nnew tolerance on the pair (1 s, 2 s) of that column: ",
    format(1e-8 * 2 + 1e-11 * max(d$delay)), " s, against a spread of ",
    "2 s\n", sep = "")

# The two rules side by side on that column, in one process, so the
# flip is a measurement rather than a citation. The old rule is written
# out here because the build that carried it is gone.
gi <- match(d$cond, sort(unique(d$cond)))
first <- match(seq_along(sort(unique(d$cond))), gi)
ref <- first[gi]
old_tol <- 1e-8 * max(abs(d$delay))
old_bad <- sort(unique(gi[abs(d$delay - d$delay[ref]) > old_tol]))
new_bad <- frmtmb.eam:::gd_varying_groups(d$delay, gi, first)
cat("\nsame column, both rules, one process:\n")
cat("  old, 1e-8 * colmax          tolerance", format(old_tol),
    " conditions flagged:", length(old_bad), "\n")
cat("  new, per pair plus a floor  tolerance",
    format(1e-8 * 2 + 1e-11 * max(abs(d$delay))),
    " conditions flagged:", length(new_bad), "\n")
cat("  centered column, new rule   conditions flagged:",
    length(frmtmb.eam:::gd_varying_groups(d2$delay, gi, first)), "\n")
cat("  centered column, old rule   conditions flagged:",
    length(sort(unique(gi[abs(d2$delay - d2$delay[ref]) >
                            1e-8 * max(abs(d2$delay))]))), "\n")
