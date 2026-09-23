# Lane wt-reunc, punch round 1 (M2): the caller's random numbers.
#
# predict_b_drawer() takes one seed per replicate. Round 1 took them
# BEFORE predict_simulate() captured the .Random.seed it restores on
# exit, so a conditional predict() moved the caller's stream by ndraws
# uniforms while returning identical predictions. A runif() after such a
# call then generated a different data set than 0.61.0 would have.
#
# The defect cannot be seen on the base build rellib-r3, which has no
# drawer at all: the failing build is round 1 of this lane and it no
# longer exists. So the round-1 ordering is RECONSTRUCTED here on the
# fixed build, by rewriting predict_simulate()'s body to move the
# capture back after the drawer and installing it in the namespace.
# That runs the defect and the fix in one process, against the base
# build's own stream position as the reference.
#
#   Rscript dev/reunc-stream.R <lane-lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(7)
G <- 10; m <- 6
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(G * m, 0, 0.8)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
nd <- d[c(1, 7, 13), c("x", "g")]

# the stream position a call leaves behind, and the next uniform after
# it, which is what a caller actually notices
after <- function(expr) {
  set.seed(4242)
  force(expr)
  list(seed = get(".Random.seed", envir = globalenv()), u = runif(1))
}

na <- after(predict(fit, newdata = nd, re_formula = NA, ndraws = 25))
cond <- after(predict(fit, newdata = nd, ndraws = 25))
set.seed(99)
pv <- predict(fit, newdata = nd, ndraws = 25)

cat("\nFIXED build\n")
cat("  stream after conditional predict() == after re_formula = NA: ",
    identical(na$seed, cond$seed), "\n", sep = "")
cat("  next runif() equal: ", identical(na$u, cond$u),
    sprintf("  (%.10f vs %.10f)\n", na$u, cond$u), sep = "")

# Reconstruct round 1: the capture and its on.exit move to AFTER the
# drawer, which is the only change. Everything else is the shipped body.
stmts <- as.list(body(frmtmb:::predict_simulate))
one <- function(pat) {
  hit <- which(vapply(stmts, function(s) {
    any(grepl(pat, paste(deparse(s), collapse = " ")))
  }, logical(1)))
  stopifnot(length(hit) == 1L)
  hit
}
i_save <- one("^saved <- get\\(\"\\.Random\\.seed\"")
i_exit <- one("^on\\.exit\\(assign\\(\"\\.Random\\.seed\", saved")
i_draw <- one("^bdraw <- predict_b_drawer\\(")
stopifnot(i_save < i_draw, i_exit < i_draw)
moved <- append(stmts[-c(i_save, i_exit)], stmts[c(i_save, i_exit)],
                after = i_draw - 2L)
f_pre <- frmtmb:::predict_simulate
body(f_pre) <- as.call(moved)
environment(f_pre) <- asNamespace("frmtmb")
assignInNamespace("predict_simulate", f_pre, ns = "frmtmb")

na1 <- after(predict(fit, newdata = nd, re_formula = NA, ndraws = 25))
cond1 <- after(predict(fit, newdata = nd, ndraws = 25))
set.seed(99)
pv1 <- predict(fit, newdata = nd, ndraws = 25)

cat("\nROUND-1 ordering, reconstructed\n")
cat("  stream after conditional predict() == after re_formula = NA: ",
    identical(na1$seed, cond1$seed), "\n", sep = "")
cat("  next runif() equal: ", identical(na1$u, cond1$u),
    sprintf("  (%.10f vs %.10f)\n", na1$u, cond1$u), sep = "")
cat("  predictions identical to the fixed build: ",
    identical(pv, pv1), "\n", sep = "")
cat("  uniforms consumed by the drawer (ndraws = 25): ",
    sum(na1$seed != cond1$seed) > 0, "\n", sep = "")
