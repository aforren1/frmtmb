# What the fit-time standard-error check costs: frm() with the check
# (default) against frm(control = frmtmb_control(check_se = "ignore")),
# interleaved in one process, each arm grown past 1.2 s per round, the
# minimum of 5 rounds. The control arm is the same build. A third arm
# repeats the default, so the ratio it reports is the noise floor.
#   Rscript dev/nanse-cost2.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
set.seed(1)
sleep <- lme4::sleepstudy
epi <- brms::epilepsy
n <- 2000
dg <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:50, n, TRUE)))
dg$y <- rpois(n, exp(0.2 + 0.3 * dg$x + rnorm(50, 0, 0.5)[dg$g]))
dl <- data.frame(x = runif(300, 0, 3))
dl$y <- 2 * exp(-0.3 * dl$x) + rnorm(300, 0, 0.2)
models <- list(
  lm = function(ctl) frm(Reaction ~ Days, data = sleep, control = ctl),
  lmm_slope = function(ctl) frm(Reaction ~ Days + (Days | Subject),
                                data = sleep, control = ctl),
  pois_glmm = function(ctl) frm(count ~ zAge + zBase * Trt + (1 | patient),
                                data = epi, family = poisson(),
                                control = ctl),
  pois_glmm_2000 = function(ctl) frm(y ~ x + z + (1 | g), data = dg,
                                     family = poisson(), control = ctl),
  dist_sigma = function(ctl) frm(bf(Reaction ~ Days, sigma ~ Days),
                                 data = sleep, control = ctl),
  nl = function(ctl) frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                         data = dl, start = list(beta = c(2, -0.3)),
                         control = ctl),
  negbin = function(ctl) frm(count ~ zAge + zBase * Trt + (1 | patient),
                             data = epi, family = negbinomial(),
                             control = ctl)
)
on <- frmtmb_control()
off <- frmtmb_control(check_se = "ignore")
per_call <- function(f, ctl) {
  reps <- 0L
  t0 <- proc.time()[[3]]
  while (proc.time()[[3]] - t0 < 1.2) {
    suppressWarnings(f(ctl))
    reps <- reps + 1L
  }
  (proc.time()[[3]] - t0) / reps
}
for (nm in names(models)) {
  a <- b <- cc <- numeric(5)
  for (r in 1:5) {
    a[r] <- per_call(models[[nm]], on)
    b[r] <- per_call(models[[nm]], off)
    cc[r] <- per_call(models[[nm]], on)
  }
  cat(sprintf(paste0("%-15s check %.4fs  ignore %.4fs  ratio %.3f  ",
                     "control ratio (check vs check) %.3f\n"),
              nm, min(a), min(b), min(a) / min(b), min(cc) / min(a)))
}
