# Lane slopetrig, item 5033: what the ratio trigger costs. For each
# design, the time of the default fit and of autoscale = FALSE, each a
# block of fits grown past 1.2 s, the minimum of 5 rounds with the arms
# interleaved (proc.time() ticks at 10 ms), and a CONTROL arm that runs
# autoscale = FALSE again, whose ratio to the first must be near 1.
# Also the shortfall of each arm against autoscale = TRUE and whether
# the default engaged.
# Usage: Rscript slopetrig-cost.R base|lane
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(lme4)
})
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "\n")
unc <- function(seed, c0, s, n_g = 20, n_per = 15, fam = "gau") {
  set.seed(seed)
  n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  z <- rnorm(n)
  eta <- 0.5 + rnorm(n_g, 0, 0.5)[g] + rnorm(n_g, 0, 0.5)[g] * z + 0.3 * z
  y <- if (fam == "gau") eta + rnorm(n) else rpois(n, exp(eta))
  data.frame(g = g, x = c0 + s * z, y = y)
}
binlev <- function(seed, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = 300)
  z <- rnorm(300)
  x <- c0 + trt + 0.01 * z
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + (0.3 + rnorm(20, 0, 0.5)[g]) * (x - c0)
  data.frame(g = g, x = x, y = eta + rnorm(300))
}
ss <- lme4::sleepstudy
designs <- list(
  list(id = "sleepstudy", d = ss, fo = Reaction ~ Days + (Days | Subject),
       fam = gaussian()),
  list(id = "unc c0=1 ratio3", d = unc(1, 1, 0.354), fo = y ~ x + (1 + x | g),
       fam = gaussian()),
  list(id = "unc c0=10 ratio2.2", d = unc(2, 10, 5), fo = y ~ x + (1 + x | g),
       fam = gaussian()),
  list(id = "binlev c0=1", d = binlev(1, 1), fo = y ~ x + (1 + x | g),
       fam = gaussian()),
  list(id = "age 40+10z", d = unc(1, 40, 10), fo = y ~ x + (1 + x | g),
       fam = gaussian()),
  list(id = "year seed2", d = unc(2, 2000, 3), fo = y ~ 1 + (1 + x | g),
       fam = gaussian()),
  list(id = "unc c0=1000 ratio4", d = unc(3, 1000, 258), fo = y ~ 1 + (1 + x | g),
       fam = gaussian()),
  list(id = "pois 5000 x=10+3z", d = unc(1, 10, 3, 100, 50, "pois"),
       fo = y ~ x + (1 + x | g), fam = poisson()),
  list(id = "gau 20000 x=10+3z", d = unc(1, 10, 3, 400, 50),
       fo = y ~ x + (1 + x | g), fam = gaussian()))
blk <- function(f) {
  n <- 1L
  repeat {
    t <- system.time(for (i in seq_len(n)) f())[["elapsed"]]
    if (t > 1.2 || n >= 4096L) return(n)
    n <- n * 2L
  }
}
for (ds in designs) {
  ctl_f <- frmtmb_control(autoscale = FALSE)
  fD <- function() suppressWarnings(suppressMessages(
    frm(ds$fo, family = ds$fam, data = ds$d)))
  fF <- function() suppressWarnings(suppressMessages(
    frm(ds$fo, family = ds$fam, data = ds$d, control = ctl_f)))
  fT <- function() suppressWarnings(suppressMessages(
    frm(ds$fo, family = ds$fam, data = ds$d,
        control = frmtmb_control(autoscale = TRUE))))
  a <- fD()
  b <- fF()
  tt <- fT()
  n <- max(blk(fD), blk(fF))
  rD <- rF <- rC <- numeric(5)
  for (k in 1:5) {
    rD[k] <- system.time(for (i in seq_len(n)) fD())[["elapsed"]] / n
    rF[k] <- system.time(for (i in seq_len(n)) fF())[["elapsed"]] / n
    rC[k] <- system.time(for (i in seq_len(n)) fF())[["elapsed"]] / n
  }
  lt <- as.numeric(logLik(tt))
  cat(sprintf(paste0("%-20s engaged %-5s n %4d default %.4f s FALSE %.4f s",
                     " default/FALSE %.2f control %.2f ; short default %.2g",
                     " FALSE %.2g ; identical par %s\n"),
              ds$id, isTRUE(a$autoscaled), n, min(rD), min(rF),
              min(rD) / min(rF), min(rC) / min(rF),
              lt - as.numeric(logLik(a)), lt - as.numeric(logLik(b)),
              identical(a$opt$par, b$opt$par)))
}
