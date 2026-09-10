# rev-gddm: is the band under the tolerance reachable by a covariate a
# field writes, and what does the fit report when it is?
#
# The tolerance is 1e-8 times the column's own largest FINITE entry, so
# it is set by the LARGEST value in the column and applied to every row
# of it. A column whose dynamic range is wider than 1e8 therefore
# carries a tolerance larger than the small entries themselves, and a
# within-condition difference the user regards as the manipulation is
# invisible to the guard.
#
# Delays from one second to ten years are the shape used in
# intertemporal choice, where a drift-diffusion analysis of the choice
# RT is routine. 3.15e8 s is ten years; 1 s, 2 s and 3 s are the short
# arm.
#
# Seed 909, 120 rows, grid dt = 0.05, ny = 51, t_max = 2.
# Arm from GDDM_LIB; default is this review's install of the worktree.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(909)
n <- 120L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))

# Four blocks of 30 trials. The index the user built names the BLOCK,
# which is the defect this guard exists to catch: `delay` is on the
# right-hand side and is not in the index.
d$blk <- rep(1:4, each = 30L)
d$cond <- d$blk
# Blocks 1 and 2 are the short arm: delays of 1, 2 and 3 seconds, drawn
# so they are not aligned with anything. Blocks 3 and 4 are the long
# arm at ten years, which is what makes the column's maximum large.
short <- c(1, 2, 3)
d$delay <- ifelse(d$blk <= 2L, sample(short, n, replace = TRUE), 3.15e8)
tol <- 1e-8 * max(abs(d$delay))
cat("delay column: max", format(max(d$delay), scientific = TRUE),
    " tolerance", format(tol, digits = 4), "s\n")
cat("within-condition spread, per block:",
    paste(tapply(d$delay, d$blk, function(z) max(z) - min(z)),
          collapse = " "), "\n")
cat("distinct delays in block 1:",
    paste(sort(unique(d$delay[d$blk == 1L])), collapse = " "), "\n\n")

form <- bf(rt | vint(upper, cond) ~ delay, bs ~ 1, ndt ~ 1, bias = 0.5)

cat("== 1. does the guard fire?\n")
m <- tryCatch({
  frm(form, family = gddm(control = ctl), data = d, dry_run = "frame")
  ""
}, error = function(e) conditionMessage(e))
cat("  ", if (nzchar(m)) paste("REFUSED:", substr(m, 1, 80)) else
       "ACCEPTED, no error and no warning", "\n", sep = "")

cat("\n== 2. does the density read the rows the guard let through?\n")
# Two data sets the guard accepts, differing only in the delay of rows
# that are not their condition's first row. If the objective is bitwise
# equal at the SAME parameter vector, those rows never reached it.
first <- vapply(split(seq_len(n), d$cond), min, integer(1))
d2 <- d
j <- setdiff(seq_len(n), first)
d2$delay[j] <- ifelse(d2$blk[j] <= 2L, 1, 3.15e8)   # every short row -> 1 s
cat("  data set B: every non-first short-arm row set to 1 s\n")
mb <- tryCatch({
  frm(form, family = gddm(control = ctl), data = d2, dry_run = "frame")
  ""
}, error = function(e) conditionMessage(e))
cat("  guard on data set B: ",
    if (nzchar(mb)) "REFUSED" else "ACCEPTED", "\n", sep = "")

obj_of <- function(dat) frm(form, family = gddm(control = ctl),
                            data = dat, dry_run = "objective")$obj
oa <- obj_of(d)
p <- oa$par
p[] <- c(0.7, -0.3, 0.9, -1.2)[seq_along(p)]
fa <- as.numeric(oa$fn(p))
fb <- as.numeric(obj_of(d2)$fn(p))
cat(sprintf("  fn(A) %.17g\n  fn(B) %.17g\n  identical: %s\n",
            fa, fb, identical(fa, fb)))
# the control: move the FIRST row of each condition instead
d3 <- d
d3$delay[first] <- d3$delay[first] + 1
fc <- as.numeric(obj_of(d3)$fn(p))
cat(sprintf("  fn, first rows only + 1 s: %.17g   moved: %s\n",
            fc, !identical(fa, fc)))

cat("\n== 3. what the fit reports\n")
fit <- try(frm(form, family = gddm(control = ctl), data = d),
           silent = TRUE)
if (inherits(fit, "try-error")) {
  cat("  fit refused: ",
      conditionMessage(attr(fit, "condition")), "\n", sep = "")
} else {
  cf <- unlist(fixef(fit))
  cat("  coefficients:",
      paste(names(cf), sprintf("%.6g", cf), sep = "=", collapse = "  "),
      "\n")
  cat("  the delay coefficient is fitted from", length(first),
      "rows of", n, "\n")
}

cat("\n== 4. the same column after the user centers it\n")
d4 <- d; d4$delayc <- d4$delay - mean(d4$delay)
cat("  centered column max", format(max(abs(d4$delayc)), digits = 4),
    " tolerance", format(1e-8 * max(abs(d4$delayc)), digits = 4), "\n")
m4 <- tryCatch({
  frm(bf(rt | vint(upper, cond) ~ delayc, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = gddm(control = ctl), data = d4, dry_run = "frame")
  ""
}, error = function(e) conditionMessage(e))
cat("  ", if (nzchar(m4)) "REFUSED" else "ACCEPTED", "\n", sep = "")

cat("\n== 5. how tight could the constant be?\n")
# the noise the tolerance exists for, relative to the column maximum,
# over the term that produces it
for (deg in 1:4) {
  x <- rep(c(0, 0.128, 0.256, 0.512), each = 30L)
  pm <- poly(x, deg)
  gi <- match(x, sort(unique(x)))
  f1 <- match(seq_len(4L), gi)
  rel <- max(abs(pm - pm[f1[gi], , drop = FALSE]) /
               rep(apply(abs(pm), 2L, max), each = nrow(pm)))
  cat(sprintf("  poly(x, %d): worst row-vs-first-row difference is %.3e",
              deg, rel))
  cat(" of its column maximum\n")
}
cat("  the shipped constant is 1.0e-08 of the column maximum\n")
