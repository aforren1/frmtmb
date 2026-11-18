# Short check of lane squant fix round 1, item 4 follow-up: the us and cs
# fits of squant-short-compat.R stopped at code 1 with the correlation at
# its limit (cs: false convergence and a gradient warning). Is that the
# family or the data? (a) the same data and formulas under gaussian() and
# student(); (b) quantile_elf() on data with a well-identified slope sd,
# 5 seeds.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-compat2.R > dev/reviews/squant-short-compat2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  ws <- character()
  ms <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 110))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, substr(conditionMessage(m), 1, 110))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%s: ERROR %s\n", lab, substr(conditionMessage(r), 1, 300)))
    return(invisible(NULL))
  }
  vc <- tryCatch(VarCorr(r), error = function(e) NULL)
  cat(sprintf("%s: conv %s | logLik %.4f | warnings %d | messages %d | class %s\n",
              lab, r$opt$convergence, as.numeric(logLik(r)), length(ws),
              length(ms), paste(class(r), collapse = "/")))
  if (length(ws)) cat("   W:", paste(unique(ws), collapse = "\n   W: "), "\n")
  if (length(ms)) cat("   M:", paste(unique(ms), collapse = "\n   M: "), "\n")
  invisible(r)
}
# (a) the data of squant-short-compat.R
set.seed(20261117)
n <- 400
d <- data.frame(x = runif(n), x2 = runif(n),
                g = factor(sample(sprintf("g%02d", 1:25), n, TRUE)))
re0 <- rnorm(25, 0, 0.6)
re1 <- rnorm(25, 0, 0.5)
gi <- as.integer(d$g)
d$y <- 1 + 2 * d$x + re0[gi] + re1[gi] * d$x + (0.5 + d$x) * rnorm(n)
cat("== (a) same data, other families ==\n")
for (fm in c("us", "cs")) {
  f <- if (fm == "us") y ~ x + (1 + x | g) else y ~ x + cs(1 + x | g)
  run(paste(fm, "gaussian"), frm(f, family = gaussian(), data = d))
  run(paste(fm, "student"), frm(f, family = student(), data = d))
  run(paste(fm, "quantile_elf smoothing 1"),
      frm(bf(f, quantile = 0.5), family = quantile_elf(smoothing = 1),
          data = d))
}
cat("\n== (b) quantile_elf(smoothing = 1), slope sd 1.5, 40 groups of 10 ==\n")
for (s in 1:5) {
  set.seed(20261200 + s)
  ng <- 40
  d <- data.frame(x = rnorm(ng * 10), g = factor(rep(seq_len(ng), each = 10)))
  b0 <- rnorm(ng, 0, 1)
  b1 <- rnorm(ng, 0, 1.5)
  d$y <- 1 + 2 * d$x + b0[d$g] + b1[d$g] * d$x + rnorm(ng * 10)
  for (fm in c("us", "cs")) {
    f <- if (fm == "us") y ~ x + (1 + x | g) else y ~ x + cs(1 + x | g)
    run(sprintf("seed %d %s q 0.5", s, fm),
        frm(bf(f, quantile = 0.5), family = quantile_elf(smoothing = 1),
            data = d))
    run(sprintf("seed %d %s q 0.9", s, fm),
        frm(bf(f, quantile = 0.9), family = quantile_elf(smoothing = 1),
            data = d))
  }
}
cat("DONE\n")
