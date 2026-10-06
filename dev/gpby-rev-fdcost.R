# Reviewer copy of dev/gpby-fdcost.R writing its own rds files.
# The cost of fitted()'s finite-difference Est.Error, counted in MODEL
# EVALUATIONS (calls of the fitted value), which is load-independent,
# with the wall clock beside it. The cells and their seeds are lane
# wt-resmooth's (dev/resmooth-batchcost.R), so the counts compare with
# the ones it recorded. The clock instrument is the reviewer's: each
# block of repeats grows past 1.2 s, the minimum of three blocks is
# reported over the repeat count, and a CONTROL with no smooth carries
# the noise.
# Usage: Rscript dev/gpby-fdcost.R base|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")

block <- function(f, target = 1.2, rounds = 3) {
  reps <- 1L
  repeat {
    t0 <- proc.time()
    for (i in seq_len(reps)) f()
    e <- (proc.time() - t0)[["elapsed"]]
    if (e >= target) break
    reps <- reps * 2L
    if (reps > 4096L) break
  }
  best <- Inf
  for (r in seq_len(rounds)) {
    t0 <- proc.time()
    for (i in seq_len(reps)) f()
    best <- min(best, (proc.time() - t0)[["elapsed"]])
  }
  best / reps
}

# every call of the fitted value, the estimate's own included, which is
# the same one call on both arms
ns <- asNamespace("frmtmb")
cnt <- 0L
trace("fitted_point", tracer = quote(cnt <<- cnt + 1L), where = ns,
      print = FALSE)
report <- function(lab, fit, newdata, re_formula) {
  cnt <<- 0L
  fv <- fitted(fit, newdata = newdata, re_formula = re_formula)
  k <- cnt
  s <- block(function() fitted(fit, newdata = newdata,
                              re_formula = re_formula))
  cat(sprintf("%-44s evals %4d | %.4f s | Est.Error[1, 2, 1] %.10g\n",
              lab, k, s, fv[1L, 2L, 1L]))
  saveRDS(fv, sprintf("dev/gpby-rev-fdcost-%s-%s.rds", arm,
                      gsub("[^A-Za-z0-9]", "", substr(lab, 1, 2))))
}

set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
report("A gp(x) n_b=160, newdata 3 NEW positions, NA", fD,
       data.frame(x = c(-1.5, 0, 1.5)), NA)
report("B gp(x) n_b=160, in sample, NA", fD, NULL, NA)
report("C gp(x) n_b=160, newdata 3 observed, NA", fD,
       dD[c(1, 5, 9), "x", drop = FALSE], NA)
set.seed(31)
ng <- 40
per <- 20
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
report("D s(x,k=8)+(1|g) 40 lv, in sample, NULL", fB, NULL, NULL)
report("E s(x,k=8)+(1|g) 40 lv, in sample, NA", fB, NULL, NA)
set.seed(43)
fG <- suppressWarnings(frm(bf(y ~ gp(x, k = 12, c = 5 / 4)),
                          family = cumulative(), data = dD))
report("F hsgp(x) k=12, in sample, NA", fG, NULL, NA)
set.seed(59)
dH <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(ng, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
report("CONTROL x+(1|g) 40 lv, in sample, NULL", fH, NULL, NULL)
cat("DONE\n")
