# Does this arm's tmbstan sample the model? Usage:
#   Rscript tmbstan121-03-probe.R <arm>
# arm: A, Bsrc, Bbin, BbinA, C120, C120A. BbinA is the CRAN binary
# (built against StanHeaders 2.39) run with the 2.32.10 pin on the path.
#
# Three instruments, from cheapest to most direct:
#   1. the static detector frmtmb.sample ships, and its refusal;
#   2. the gradient HMC uses, rstan::grad_log_prob(), against the
#      objective's own gradient at fixed points, chain-free;
#   3. chains on two models with a known answer: an RTMB product of
#      independent normals whose posterior is exact, and kaskr/tmbstan#33's
#      own reproduction, TMB's "simple" example.
ROOT <- "C:/Users/adf44/source/r/tmbstan121-lib"
PIN <- "C:/Users/adf44/source/r/pinlib"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
arm <- commandArgs(trailingOnly = TRUE)[1]
tlib <- file.path(ROOT, if (arm == "BbinA") "Bbin" else arm)
com <- file.path(ROOT, "common")
paths <- if (arm %in% c("A", "C120A", "BbinA")) {
  c(tlib, com, PIN, USER)
} else {
  c(tlib, com, USER)
}
.libPaths(paths)

cat("ARM", arm, "\n")
pk <- c("tmbstan", "rstan", "StanHeaders", "RcppParallel", "BH", "TMB",
        "RTMB", "RcppEigen", "Rcpp", "frmtmb", "frmtmb.sample")
for (p in pk) {
  cat(sprintf("VERSION %-14s %-10s %s\n", p, format(packageVersion(p)),
              dirname(find.package(p))))
}

suppressMessages({
  library(tmbstan)
  library(RTMB)
})

# 1. The shipped detector.
broken <- frmtmb.sample:::tmbstan_build_broken()
refusal <- tryCatch({
  frmtmb.sample:::check_tmbstan_build("probe")
  "no refusal"
}, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 60)))
cat("DETECTOR tmbstan_build_broken() =", broken, "\n")
cat("DETECTOR check_tmbstan_build():", refusal, "\n")

# 2 and 3a. Independent normals with an exact posterior. Distinct means
# and scales, so a standard normal cannot match any coordinate by luck.
mu <- c(3, -2, 10, 0.5)
sdv <- c(0.5, 2, 0.1, 3)
f <- function(p) {
  -sum(dnorm(p$th, mu, sdv, log = TRUE))
}
obj <- MakeADFun(f, list(th = rep(0, 4)), silent = TRUE)

sf <- tmbstan(obj, chains = 1, iter = 20, seed = 1, refresh = 0)
set.seed(4030)
pts <- list(origin = rep(0, 4), mode = mu, mode_plus = mu + 0.3,
            random = rnorm(4, 0, 2))
for (nm in names(pts)) {
  u <- pts[[nm]]
  g_stan <- rstan::grad_log_prob(sf, u)
  g_obj <- -as.numeric(obj$gr(u))
  gap <- max(abs(as.numeric(g_stan) - g_obj)) / max(abs(g_obj), 1)
  lp_stan <- as.numeric(attr(g_stan, "log_prob"))
  lp_obj <- -obj$fn(u)
  # The defect's closed form: a standard normal kernel, gradient -u.
  gap_def <- max(abs(-u - g_obj)) / max(abs(g_obj), 1)
  cat(sprintf(
    "GRAD point=%-9s rel_gap=%.3e identical=%s defect_form_gap=%.3e lp_stan=%.6f lp_obj=%.6f\n",
    nm, gap, identical(as.numeric(g_stan), g_obj), gap_def, lp_stan,
    lp_obj))
}

fit <- tmbstan(obj, chains = 4, iter = 4000, seed = 20260917,
               refresh = 0)
dr <- as.matrix(fit)[, 1:4]
m <- colMeans(dr)
s <- apply(dr, 2, sd)
# Monte Carlo standard error of the mean from rstan's own ess.
su <- rstan::summary(fit)$summary[1:4, ]
for (i in 1:4) {
  cat(sprintf(
    "NORMAL th[%d] truth mean=%6.2f sd=%5.2f | draws mean=%9.4f sd=%7.4f | z=(mean-truth)/se_mean=%7.2f sd/truth=%.4f n_eff=%.0f\n",
    i, mu[i], sdv[i], m[i], s[i], (m[i] - mu[i]) / su[i, "se_mean"],
    s[i] / sdv[i], su[i, "n_eff"]))
}

# 3b. kaskr/tmbstan#33's reproduction. The example is copied out of the
# TMB installation, because runExample() compiles in its folder and
# that folder is in the shared user library.
ex <- file.path(ROOT, "tmbex", arm)
dir.create(ex, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(system.file("examples", package = "TMB"),
                    c("simple.cpp", "simple.R")), ex, overwrite = TRUE)
old <- setwd(ex)
suppressMessages(TMB::compile("simple.cpp"))
env <- new.env()
suppressMessages(sys.source("simple.R", envir = env))
setwd(old)
tobj <- env$obj
topt <- env$opt
rep <- TMB::sdreport(tobj)
est <- summary(rep, "fixed")
tfit <- tmbstan(tobj, chains = 1, iter = 4000, seed = 1, refresh = 0)
ts <- rstan::summary(tfit)$summary
for (nm in c("beta[1]", "beta[2]", "logsdu", "logsd0")) {
  k <- match(sub("\\[.*", "", nm), rownames(est))
  if (nm == "beta[2]") k <- k + 1
  cat(sprintf(
    "SIMPLE %-8s ML est=%9.4f se=%7.4f | draws mean=%9.4f sd=%7.4f n_eff=%.0f\n",
    nm, est[k, 1], est[k, 2], ts[nm, "mean"], ts[nm, "sd"],
    ts[nm, "n_eff"]))
}
cat("DONE", arm, "\n")
