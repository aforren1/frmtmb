# Lane optrap, punch round 2: what predictor_units() costs and gives on
# designs far above the least-squares cap of 5e6 cells, where it fits
# evenly spaced rows up to the cap: n = 1e6 by 50 dense columns, and
# n = 1e5 with a 500-level factor (a 1e5 by 500 design, mostly zeros).
# y = 1 + 1e4 z1 + b_g + e, residual spread 1. Frame only.
#   Rscript dev/optrap-p2-lscost.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pu <- frmtmb:::predictor_units
cat("lib", dirname(find.package("frmtmb")), "\n")
one <- function(lab, d, fo) {
  fr <- frm(fo, data = d, dry_run = "frame")
  X <- fr$linpreds[[1]]$X
  t <- system.time(u <- pu(fr))[["elapsed"]]
  cat(sprintf("%-22s X %s %d x %d | unit %s | predictor_units() %.2f s\n",
              lab, class(X)[1], nrow(X), ncol(X),
              if (length(u)) paste(u, collapse = ",") else "none", t))
}
set.seed(1)
n <- 1e6
Z <- matrix(rnorm(n * 49), n, 49)
colnames(Z) <- paste0("z", 1:49)
d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
one("1e6 x 50 dense", d, as.formula(paste("y ~", paste(colnames(Z),
                                                     collapse = " + "),
                                         "+ (1 | g)")))
rm(Z, d)
gc(verbose = FALSE)
set.seed(2)
n <- 1e5
d <- data.frame(z1 = rnorm(n), f = factor(sample(500, n, TRUE)),
                g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(500)[d$f] + rnorm(1000, 0, 0.5)[d$g] +
  rnorm(n)
one("1e5 x 500 factor", d, y ~ z1 + f + (1 | g))
