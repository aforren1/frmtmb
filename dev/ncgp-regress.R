# Non-gp fits and the approximate gp(k =), digit for digit against rellib-r16:
# same tape, same numbers. Each line holds the tape's node count and %.17g of the
# objective, the outer parameters and a sum over fitted(). Seeds in place. Run
# pinned to one core type (rtmb-pitfalls item 21):
#   Rscript dev/ncgp-regress.R <lib or base> > dev/ncgp-out/regress-<arm>.txt
args <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (args[1] != "base") libs <- c(args[1], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261107)
n <- 200L
d <- data.frame(x = runif(n, 0, 6), z = rnorm(n),
                g = factor(rep(1:10, each = 20L)),
                t = rep(1:20, 10L))
d$y <- sin(d$x) + 0.5 * d$z + rnorm(10)[d$g] * 0.5 + rnorm(n, 0, 0.4)
d$yu <- d$y
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$z))
d$y2 <- 0.5 * d$y + rnorm(n, 0, 0.5)
d$bin <- rbinom(n, 1, plogis(0.3 * d$z))
d$xs <- sort(d$x)
sp <- list(
  lm = list(y ~ x + z, gaussian()),
  re = list(y ~ z + (1 | g), gaussian()),
  res = list(y ~ z + (1 + z | g), gaussian()),
  smooth = list(y ~ s(x, k = 6), gaussian()),
  sigma = list(bf(y ~ z, sigma ~ x), gaussian()),
  pois = list(cnt ~ z + (1 | g), poisson()),
  bern = list(bin ~ z + (1 | g), bernoulli()),
  ar = list(y ~ z + ar(time = t, gr = g), gaussian()),
  mv = list(mvbf(bf(y ~ z), bf(y2 ~ z), rescor = TRUE), NULL),
  mix = list(bf(y ~ 1, mu1 ~ z, mu2 ~ 1), mixture(gaussian, gaussian)),
  hsgp_unsorted = list(y ~ gp(x, k = 12), gaussian()),
  hsgp_sorted = list(bf(y ~ gp(xs, k = 12)), gaussian()),
  hsgp_by = list(y ~ gp(x, k = 8, by = g), gaussian()),
  hsgp_2d = list(y ~ gp(x, z, k = 5), gaussian()),
  hsgp_pois = list(cnt ~ gp(x, k = 8), poisson())
)
for (nm in names(sp)) {
  s <- sp[[nm]]
  f <- if (is.null(s[[2]])) {
    suppressWarnings(frm(s[[1]], data = d))
  } else {
    suppressWarnings(frm(s[[1]], data = d, family = s[[2]]))
  }
  tp <- tryCatch(nrow(RTMB::GetTape(f$obj, "ADFun")$data.frame()), error = function(e) NA)
  fv <- sum(as.matrix(suppressWarnings(fitted(f)))[, 1L])
  cat(sprintf("%s nodes %s obj %.17g fitted %.17g par %s\n", nm, tp, f$opt$objective, fv,
              paste(sprintf("%.17g", f$opt$par), collapse = " ")))
}
