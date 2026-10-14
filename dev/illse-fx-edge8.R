# Lane illse fix round 1: degree-8 fits with (1 | g) and the sd at zero
# (dev/illse-rev-edge4.R design): which SEs are kept, and how they
# compare with lm()'s (ML scale), which is the conditional covariance
# with the sd at its boundary.
#   Rscript dev/illse-fx-edge8.R base|lane|dev seed...
args <- commandArgs(TRUE)
arm <- args[1]
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
deg <- 8
for (seed in as.integer(args[-1])) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  f <- suppressWarnings(frm(fo, data = d))
  lmf <- stats::lm(stats::as.formula(paste("y ~", paste(xn, collapse = "+"))),
                   data = d)
  se_l <- sqrt(diag(vcov(lmf)) * (n - deg - 1) / n)
  se_f <- suppressWarnings(sqrt(diag(vcov(f))))
  l <- ns$sdr_of(f)$se_lost
  cat(sprintf("seed %d | lost %s\n", seed,
              paste(names(l), l, collapse = ", ")))
  cat("  SE / lm:", formatC(se_f[seq_along(se_l)] / se_l, digits = 6,
                            format = "g"), "\n")
}
