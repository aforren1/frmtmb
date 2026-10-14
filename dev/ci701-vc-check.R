# Lane ci701: confint_varcorr() against VarCorr() on the y ~ x +
# (1 + x | g) fits of dev/ci701-vc-spread.R whose correlation the
# standard-error check removed, and how stable that outcome is when x is
# moved by a few ulps (seed 3000 + k).
#
#   Rscript dev/ci701-vc-check.R <lib|base> <seeds, comma separated> <K>
args <- commandArgs(TRUE)
lib <- args[1]
seeds <- as.integer(strsplit(args[2], ",")[[1]])
K <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
for (s in seeds) {
  set.seed(s)
  ng <- 20
  g <- factor(rep(seq_len(ng), each = 10))
  d0 <- data.frame(g = g, x = stats::rnorm(200))
  d0$y <- 1 + 0.5 * d0$x + stats::rnorm(ng, 0, 0.7)[g] + stats::rnorm(200)
  for (k in 0:K) {
    d <- d0
    if (k > 0) {
      set.seed(3000 + k)
      d$x <- d$x * (1 + sample(0:4, nrow(d), TRUE) * 2^-52)
    }
    fit <- suppressMessages(suppressWarnings(
      frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)))
    lost <- frmtmb:::sdr_of(fit)$se_lost
    ci <- suppressWarnings(confint_varcorr(fit))
    vc <- suppressWarnings(VarCorr(fit))
    sdr <- ci[ci$type == "sd", ]
    if (k == 0 && s == seeds[1]) print(str(vc))
    cat(sprintf("RES seed=%d k=%d code=%d lost=%s sd_lwr=%s sd_upr=%s\n",
                s, k, fit$opt$convergence,
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "none",
                paste(sprintf("%.5f", sdr$lwr), collapse = ","),
                paste(sprintf("%.5f", sdr$upr), collapse = ",")))
  }
}
