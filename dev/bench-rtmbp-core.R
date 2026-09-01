## Head-to-head RTMB vs RTMBp on a hand-rolled objective.
##
## One configuration per R process (RTMB and RTMBp both define an S4 class
## "advector"; attaching both in one session makes method dispatch
## ambiguous, so they are never loaded together).
##
## usage: Rscript dev/bench-rtmbp-core.R <backend> <threads> <shape> [autopar] [mode] [pass]
##   backend  RTMB | RTMBp
##   threads  integer, passed to TMB::openmp(n, DLL = backend)
##   shape    A (InstEval LMM) | B (Poisson GLMM n = 3e5)
##   autopar  TRUE (default) | FALSE
##   mode     full (default) | calls  (calls = tape + identity + per-call
##            only, cheap enough to run in several interleaved passes so
##            thermal drift on this laptop does not alias onto one config)
##   pass     integer label, appended to the output file name
##
## Threads and autopar are set BEFORE MakeADFun: TMBad's automatic
## parallelization splits the tape at tape time, so every configuration
## re-tapes.
args <- commandArgs(trailingOnly = TRUE)
backend <- args[[1]]
threads <- as.integer(args[[2]])
shape   <- args[[3]]
autopar <- if (length(args) >= 4) as.logical(args[[4]]) else TRUE
mode    <- if (length(args) >= 5) args[[5]] else "full"
pass    <- if (length(args) >= 6) args[[6]] else ""

.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/libp",
            .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(lme4)
  library(backend, character.only = TRUE)
})
source("dev/bench-rtmbp-data.R")

set_threads <- function() {
  TMB::openmp(threads, autopar = autopar, DLL = backend)
}
cfg_before <- set_threads()
cfg <- TMB::config(DLL = backend)

tag <- sprintf("%s-t%d-%s-ap%s%s%s", backend, threads, shape, autopar,
               if (mode == "calls") "-calls" else "",
               if (nzchar(pass)) paste0("-p", pass) else "")
cat(sprintf("=== %s ===  nthreads=%s autopar=%s\n", tag,
            cfg$nthreads, cfg$autopar))

ref <- readRDS("dev/rtmbp-ref.rds")
popt <- ref[[paste0("opt_", shape)]]

dat <- build_shape(shape)
f <- make_f(shape, dat)

R <- list(tag = tag, backend = backend, threads = threads, shape = shape,
          autopar = autopar, nthreads_cfg = cfg$nthreads,
          autopar_cfg = cfg$autopar,
          Rver = R.version.string,
          pkgver = as.character(packageVersion(backend)))

mk <- function() {
  set_threads()
  MakeADFun(f, dat$parameters, random = dat$random, silent = TRUE)
}

## --- (a) tape build time ----------------------------------------------
nt <- 5L
t_tape <- numeric(nt)
for (i in seq_len(nt)) { gc(); t_tape[i] <- tm(obj <- mk()) }
R$tape <- t_tape
cat(sprintf("tape:  min %.3f med %.3f max %.3f\n",
            min(t_tape), median(t_tape), max(t_tape)))

## --- identity: fn / gr at fixed reference points ----------------------
set.seed(11)
np <- length(popt)
pts <- c(list(popt), lapply(1:3, function(i) popt + c(0.01, -0.02, 0.015,
                                                      0.005, -0.01)[seq_len(np)] * i))
obj$fn(popt)                              # warm the inner solve
R$id_fn <- vapply(pts, function(q) obj$fn(q), 0)
R$id_gr <- lapply(pts, function(q) as.numeric(obj$gr(q)))
R$id_pts <- pts

## --- (b) per-call cost at the optimum ---------------------------------
n <- 25L
set.seed(11)
jit <- lapply(seq_len(n), function(i) popt + rnorm(np, 0, 0.02))
obj$fn(popt)
t_fn <- t_grc <- numeric(n)
for (i in seq_len(n)) {
  q <- jit[[i]]
  t_fn[i]  <- tm(obj$fn(q))     # inner Newton solve at a new point
  t_grc[i] <- tm(obj$gr(q))     # adjoint, last.par == q
}
set.seed(23)
t_grf <- numeric(n)
for (i in seq_len(n)) {
  q <- popt + rnorm(np, 0, 0.02)
  t_grf[i] <- tm(obj$gr(q))     # fresh point: solve + adjoint
}
R$fn <- t_fn; R$gr_cached <- t_grc; R$gr_fresh <- t_grf

## CPU vs elapsed over a block of fn+gr calls: cpu/elapsed > 1 means
## OpenMP threads really ran, = 1 means the tape stayed serial.
R$callblock <- tmc(for (i in seq_len(n)) { q <- jit[[i]]; obj$fn(q); obj$gr(q) })
cat(sprintf("block: elapsed %.2f cpu %.2f  ratio %.2f\n",
            R$callblock[["elapsed"]], R$callblock[["cpu"]],
            R$callblock[["cpu"]] / R$callblock[["elapsed"]]))
cat(sprintf("fn:    min %.4f med %.4f max %.4f\n", min(t_fn), median(t_fn), max(t_fn)))
cat(sprintf("gr(c): min %.4f med %.4f max %.4f\n", min(t_grc), median(t_grc), max(t_grc)))
cat(sprintf("gr(f): min %.4f med %.4f max %.4f\n", min(t_grf), median(t_grf), max(t_grf)))

## --- (c) full nlminb fit ----------------------------------------------
if (mode == "full") {
nf <- 5L
t_fit <- numeric(nf); c_fit <- numeric(nf)
pars <- vector("list", nf); nlls <- numeric(nf)
its <- numeric(nf); fev <- numeric(nf); gev <- numeric(nf)
for (i in seq_len(nf)) {
  gc()
  o2 <- mk()
  z <- tmc(op <- nlminb(o2$par, o2$fn, o2$gr))
  t_fit[i] <- z[["elapsed"]]; c_fit[i] <- z[["cpu"]]
  pars[[i]] <- op$par; nlls[i] <- op$objective
  its[i] <- op$iterations
  fev[i] <- op$evaluations[["function"]]; gev[i] <- op$evaluations[["gradient"]]
  if (i == nf) obj_fit <- o2
}
R$fit <- t_fit; R$fit_cpu <- c_fit; R$fit_par <- pars[[nf]]; R$fit_nll <- nlls
R$fit_iter <- its; R$fit_fev <- fev; R$fit_gev <- gev
cat(sprintf("fit:   min %.2f med %.2f max %.2f  cpu/el %.2f  (%d it, %d fn, %d gr, nll %.6f)\n",
            min(t_fit), median(t_fit), max(t_fit),
            median(c_fit) / median(t_fit), its[nf], fev[nf], gev[nf], nlls[nf]))

## --- (d) sdreport ------------------------------------------------------
ns <- 5L
t_sd <- numeric(ns)
for (i in seq_len(ns)) {
  gc()
  t_sd[i] <- tm(sdr <- sdreport(obj_fit))
}
R$sdreport <- t_sd
ss <- summary(sdr, "fixed")
R$sd_est <- ss[, 1]; R$sd_se <- ss[, 2]
cat(sprintf("sdrep: min %.2f med %.2f max %.2f\n",
            min(t_sd), median(t_sd), max(t_sd)))
}

dir.create("dev/rtmbp-results", showWarnings = FALSE)
saveRDS(R, file.path("dev/rtmbp-results", paste0(tag, ".rds")))
cat("wrote dev/rtmbp-results/", tag, ".rds\n", sep = "")
