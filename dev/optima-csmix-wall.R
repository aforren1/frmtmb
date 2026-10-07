# Lane optima, item 3: are the cs() ordinal mixtures that stop with
# "false convergence (8)" against the crossing wall? For each fit of
# dev/optima-csmix.R's mix_cum_sratio, the smallest gap between a
# row's two thresholds in its own category in the cumulative component
# (tau_k - cs_k z - eta, on the link scale), against the optimizer code
# and max |gradient|.
#   Rscript dev/optima-csmix-wall.R base|lane [seeds]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:20
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
rows <- list()
for (s in seeds) {
  d <- mk(s, mix = TRUE)
  f <- tryCatch(suppressWarnings(frm(bf(y ~ x + cs(z)),
                                     family = mixture(cumulative(),
                                                      sratio()),
                                     data = d)),
                error = function(e) NULL)
  if (is.null(f)) {
    rows[[length(rows) + 1L]] <- data.frame(seed = s, code = NA,
                                            grad = NA, min_gap = NA)
    next
  }
  pl <- f$estimates
  tr <- pl$tau_raw1
  tau <- cumsum(c(tr[1], exp(tr[-1])))
  eta <- pl$beta[1] * d$x
  A <- outer(-eta, tau, "+") - outer(d$z, pl$bcs3)
  lo <- cbind(-Inf, A)
  up <- cbind(A, Inf)
  i <- cbind(seq_len(nrow(d)), d$y)
  gap <- up[i] - lo[i]
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, code = f$opt$convergence,
    grad = max(abs(f$obj$gr(f$opt$par))), min_gap = min(gap))
}
X <- do.call(rbind, rows)
print(X)
cat("code 1 fits:", sum(X$code %in% 1), "; of them with min gap < 1e-3:",
    sum(X$code %in% 1 & X$min_gap < 1e-3), "\n")
cat("code 0 fits:", sum(X$code %in% 0), "; min gaps:",
    paste(signif(X$min_gap[X$code %in% 0], 3), collapse = " "), "\n")
