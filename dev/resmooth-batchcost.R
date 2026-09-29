# Lane wt-resmooth, nits round. The cost of the finite-difference
# Est.Error route counted in MODEL EVALUATIONS, which is
# load-independent, beside the wall clock. The clock instrument is the
# reviewer's: proc.time() ticks at 10 ms here, so each arm grows a block
# of repeats past 1.2 s and reports the minimum of three blocks over the
# repeat count, with a CONTROL that must report a ratio near 1.
# Run it with RESMOOTH_LANE_LIB set to the library to measure, into
# dev/resmooth-batchcost-before.txt and -after.txt.
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
LANE <- if (nzchar(Sys.getenv("RESMOOTH_LANE_LIB"))) {
  Sys.getenv("RESMOOTH_LANE_LIB")
} else "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(if (!base) LANE,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")

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

# Count the model evaluations fit_fd_se() makes for the b_idx and the
# batches the shipped route would use, by handing it an `f` that counts.
# The count is exact and does not move with machine load.
nev <- function(fit, newdata, re_formula) {
  b_idx <- {
    want <- frmtmb:::smooth_b_idx(fit)
    if (frmtmb:::re_form_keeps(re_formula)) {
      want <- sort(unique(c(want, frmtmb:::re_governed_b(fit))))
    }
    if (!length(want)) NULL else {
      used <- frmtmb:::re_used_b(fit, newdata, NULL, FALSE)
      if (is.null(used)) want else intersect(want, used)
    }
  }
  bt <- if (length(b_idx)) {
    frmtmb:::re_b_batches(fit, newdata, NULL, FALSE, b_idx)
  }
  n <- 0L
  f <- function(x) {
    n <<- n + 1L
    frmtmb:::fitted_point(x, newdata, re_formula)
  }
  invisible(frmtmb:::fit_fd_se(fit, f, b_idx = b_idx, b_batch = bt))
  c(evals = n, b = length(b_idx),
    batches = if (is.null(bt)) NA_integer_ else length(bt))
}

report <- function(lab, fit, newdata, re_formula) {
  k <- nev(fit, newdata, re_formula)
  s <- block(function() fitted(fit, newdata = newdata,
                              re_formula = re_formula))
  cat(sprintf(paste0("%-46s evals %4d | b %3d | batches %3s | ",
                     "%.4f s\n"),
              lab, k[["evals"]], k[["b"]],
              if (is.na(k[["batches"]])) "NULL" else k[["batches"]], s))
}

# A: the reviewer's headline cell. An exact gp() at THREE NEW positions:
# every kriging row loads every one of the 160 columns.
set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
report("A gp(x) n_b=160, newdata 3 NEW positions, NA", fD,
       data.frame(x = c(-1.5, 0, 1.5)), NA)
# B: the same fit IN SAMPLE, where the gp design is an indicator and one
# batch is exact. This is the cell the minimal fix is for.
report("B gp(x) n_b=160, in sample, NA", fD, NULL, NA)
# C: and at three OBSERVED positions, the same indicator fast path
report("C gp(x) n_b=160, newdata 3 observed, NA", fD,
       dD[c(1, 5, 9), "x", drop = FALSE], NA)

# D, E: the reviewer's smooth cells. Every row loads every basis column,
# so one pair per column is irreducible under finite differences.
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
# F: an hsgp() fit, whose basis is dense like a smooth's
set.seed(43)
dG <- dD
fG <- suppressWarnings(frm(bf(y ~ gp(x, k = 12, c = 5 / 4)),
                          family = cumulative(), data = dG))
report("F hsgp(x) k=12, in sample, NA", fG, NULL, NA)

# CONTROL: no smooth and no gp, so nothing this round can touch it
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
