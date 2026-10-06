# Lane fixes, punch round 2, m11: the reviewer's lane-only silent case
# F11 seed 8, bf(y ~ s(x1), sigma ~ s(x2)) on gamSim eg 6 with the
# reviewer's data builder (dev/fixes-rev2-smooth.R): which optimizer
# setting reaches the better optimum the base build finds.
#   Rscript dev/fixes-p2-f11.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
gs <- function(eg, seed, n = 200) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = eg, n = n, verbose = FALSE))
  d$z <- runif(nrow(d))
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
  d$gr <- factor(rep(1:20, length.out = nrow(d)))
  d$y <- d$y + 0.6 * stats::rnorm(20)[d$gr]
  d
}
for (seed in c(8, 1, 2, 3)) {
  d <- gs(6, seed)
  fo <- bf(y ~ s(x1), sigma ~ s(x2))
  one <- function(lab, ctl = frmtmb_control()) {
    f <- suppressWarnings(frm(fo, data = d, control = ctl))
    cat(sprintf("seed %d %-34s conv %d logLik %.6f theta %s\n", seed, lab,
                f$opt$convergence, as.numeric(logLik(f)),
                paste(sprintf("%.2f", f$estimates$theta), collapse = " ")))
  }
  one("default")
  one("restarts = 3", frmtmb_control(restarts = 3))
  if (exists("smooth_fx_units", asNamespace("frmtmb"))) {
    real <- get("smooth_fx_units", asNamespace("frmtmb"))
    assignInNamespace("smooth_fx_units", function(...) NULL, ns = "frmtmb")
    one("no smooth units")
    assignInNamespace("smooth_fx_units", real, ns = "frmtmb")
  }
}
