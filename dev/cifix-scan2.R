# Scan seeds and sizes around test-difference.R:490 for a fit on which
# frm_curve(contrast = ) refuses at its covariance check, and record why.
# Usage: Rscript dev/cifix-scan2.R <lib or "base"> <seed from> <seed to>
#        <out tsv>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
seeds <- as.integer(args[2]):as.integer(args[3])
out <- args[4]
blas <- system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]]
rows <- list()
for (n in c(100L, 120L, 140L)) for (s in seeds) {
  set.seed(s)
  d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                  fac = factor(rep(c("A", "B"), length.out = n)))
  d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
    stats::rnorm(n, 0, 0.3)
  fit <- tryCatch(suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)),
        family = stats::gaussian(), data = d)), error = function(e) NULL)
  if (is.null(fit)) next
  gx <- d$x[-1] - diff(d$x) / 2
  A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
  B <- data.frame(x = gx, fac = factor("B", levels = levels(d$fac)))
  qmin <- Inf
  emax <- 0
  nlost <- 0L
  for (g in list(A, B)) {
    lb <- suppressWarnings(frm_lp_basis(fit, newdata = g))
    C <- as.matrix(lb$A)
    q <- rowSums((C %*% lb$V) * C)
    qmin <- min(qmin, q)
    emax <- max(emax, lb$extra_var)
    nlost <- nlost + sum(lb$se_nonest)
  }
  msg <- tryCatch({
    suppressWarnings(frm_curve(fit, newdata = A, contrast = B,
                               simultaneous = FALSE))
    "ok"
  }, error = function(e) gsub("[\t\n]", " ", conditionMessage(e)))
  lost <- frmtmb:::sdr_of(fit)$se_lost
  rows[[length(rows) + 1L]] <- data.frame(
    n = n, seed = s, th1 = fit$estimates$theta[1],
    th2 = fit$estimates$theta[2], qmin = qmin, emax = emax,
    nlost_rows = nlost, se_lost = paste(names(lost), lost, collapse = ";"),
    msg = substr(msg, 1, 140), blas = blas)
  cat(n, s, signif(qmin, 4), substr(msg, 1, 60), "\n")
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE,
                   row.names = FALSE)
