# Time of one-block quadrature fits, base against lane, and where the
# lane's time goes (verbose stages). Interleave arms by running this
# once per library in parallel processes; each design is fitted `reps`
# times and the minimum kept.
#   Rscript dev/quadseam-timing.R <lib|base> [reps]
# log dev/quadseam-log/timing-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
reps <- if (length(a) > 1) as.integer(a[2]) else 3L
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
designs <- list()
set.seed(1)
n <- 200
xs <- -abs(rnorm(n)) * 3
y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dd <- data.frame(y = y, xs = xs)
set.seed(99)
dd$g <- factor(rep(1:20, length.out = nrow(dd)))
designs$skewn <- list(bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1),
                      skew_normal(), dd)
set.seed(3)
d1 <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
d1$y <- rbinom(120, 1, plogis(0.2 + 0.5 * d1$x))
designs$bern_sd0 <- list(bf(y ~ x + (1 | g)), bernoulli(), d1)
set.seed(11)
d2 <- data.frame(g = factor(rep(1:40, each = 5)), x = rnorm(200))
d2$y <- rbinom(200, 1, plogis(-0.3 + 0.7 * d2$x + rnorm(40, 0, 1)[d2$g]))
designs$bern_sd1 <- list(bf(y ~ x + (1 | g)), bernoulli(), d2)
set.seed(12)
d3 <- data.frame(g = factor(rep(1:60, each = 6)), x = rnorm(360))
d3$y <- rpois(360, exp(0.2 + 0.5 * d3$x + rnorm(60, 0, 0.5)[d3$g]))
designs$pois <- list(bf(y ~ x + (1 | g)), poisson(), d3)
for (nm in names(designs)) {
  ds <- designs[[nm]]
  tt <- numeric(reps)
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    f <- suppressMessages(suppressWarnings(
      frm(ds[[1]], family = ds[[2]], data = ds[[3]], quadrature = TRUE)))
    tt[r] <- proc.time()[["elapsed"]] - t0
  }
  cat(sprintf("TIME %-9s min %.2f max %.2f logLik %.10f\n", nm, min(tt),
              max(tt), as.numeric(logLik(f))))
}
if (a[1] != "base") {
  for (nm in names(designs)) {
    ds <- designs[[nm]]
    cat("VERBOSE", nm, "\n")
    invisible(suppressWarnings(
      frm(ds[[1]], family = ds[[2]], data = ds[[3]], quadrature = TRUE,
          verbose = TRUE)))
  }
}
cat("END\n")
# tapes per fit: every RTMB::MakeADFun() call the fit makes, counted
# with trace() (the lane's settle and boundary probes each cost one)
cnt <- 0L
suppressMessages(trace("MakeADFun", where = asNamespace("RTMB"),
                       tracer = quote(cnt <<- cnt + 1L), print = FALSE))
for (nm in names(designs)) {
  ds <- designs[[nm]]
  cnt <- 0L
  invisible(suppressMessages(suppressWarnings(
    frm(ds[[1]], family = ds[[2]], data = ds[[3]], quadrature = TRUE))))
  cat(sprintf("TAPES %-9s %d\n", nm, cnt))
}
cat("END2\n")
