# Two questions the timing raised, answered in one process.
#
# 1. Does repeating `obj$gr()` at ONE parameter measure the gradient an
#    optimizer performs? `tests/testthat/helper-scale.R`'s `scale_grad()`,
#    which produced the "one gradient" column of `dev/scale-findings.md`,
#    does exactly that. TMB's Laplace objective caches the inner Newton
#    solution on the parameter, so calls two and later at one point skip
#    work an optimizer never skips. Both instruments are run here on the
#    same objects, interleaved.
#
# 2. Why is the whole fit 3.4x apart when one gradient is 1.3x apart?
#    The nlminb evaluation counts say.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB); library(frmtmb); library(frmtmb.ode)
})
args <- commandArgs(trailingOnly = TRUE)
NS <- if (length(args) >= 1L) as.integer(args[[1]]) else 25L
ROUNDS <- if (length(args) >= 2L) as.integer(args[[2]]) else 4L
# "grad" skips the fits, so that the gradient comparison can be run at
# the Phase 0 size without paying an hour for the frm_ode() fit
MODE <- if (length(args) >= 3L) args[[3]] else "all"
# a third argument of "gr" skips the fits, which at 100 subjects means
# skipping a 3716 s frm_ode() call
FITS <- !(length(args) >= 3L && identical(args[[3]], "gr"))

TR <- list(ka = 1.0, ke = 0.15, V = 20, sigma = 0.3, sd_lka = 0.3,
           sd_lke = 0.25, amt = 100, ii = 12)
ode_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
ode_doses <- data.frame(time = c(0, 12), state = "depot",
                        value = TR$amt, ii = c(12, 12),
                        addl = c(0L, 12L), ss = c(TRUE, FALSE))
conc_ss <- function(t, ka, ke, V, amt, ii) {
  u <- t %% ii
  amt * ka / (V * (ka - ke)) *
    (exp(-ke * u) / (1 - exp(-ke * ii)) -
       exp(-ka * u) / (1 - exp(-ka * ii)))
}
set.seed(20260908L)
lka <- log(TR$ka) + stats::rnorm(NS, 0, TR$sd_lka)
lke <- log(TR$ke) + stats::rnorm(NS, 0, TR$sd_lke)
tt <- 144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)
d <- data.frame(id = factor(rep(seq_len(NS), each = length(tt))),
                time = rep(tt, times = NS))
i <- as.integer(d$id)
d$conc <- conc_ss(d$time, exp(lka[i]), exp(lke[i]), TR$V, TR$amt,
                  TR$ii) + stats::rnorm(nrow(d), 0, TR$sigma)
ST <- list(beta = c(log(0.8), log(0.2), log(15)))

f_ode <- function() {
  doses <- ode_doses
  bf(conc ~ frm_ode(ode_dyn, init = list(0, 0), times = time,
                    parms = list(exp(lka), exp(lke), exp(lV)),
                    group = id, states = c("depot", "central"),
                    output = "central", events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
f_l20 <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses, n_ss = 20L),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
f_inf <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
forms <- list(ode = f_ode(), lin20 = f_l20(), linInf = f_inf())
objs <- lapply(forms, function(f)
  frm(f + gaussian(), data = d, start = ST, dry_run = "objective")$obj)
p0 <- objs$ode$par
set.seed(1L)
PJ <- matrix(rnorm(64 * length(p0), 0, 0.02), 64, byrow = TRUE) +
  rep(p0, each = 64)

el <- function(e) {
  t0 <- proc.time()[["elapsed"]]
  force(e)
  proc.time()[["elapsed"]] - t0
}
one_pt <- function(o, n) el(for (k in seq_len(n)) o$gr(p0))
many_pt <- function(o, n)
  el(for (k in seq_len(n)) o$gr(PJ[(k - 1L) %% 64L + 1L, ]))
grow <- function(o, f) {
  o$gr(p0)
  n <- 1L
  repeat {
    if (f(o, n) > 2.5 || n >= 1024L) return(n)
    n <- n * 2L
  }
}
cat("subjects:", NS, " rows:", nrow(d), " rounds:", ROUNDS, "\n\n")
res <- matrix(NA_real_, ROUNDS, 6)
colnames(res) <- c("ode.1pt", "ode.npt", "lin20.1pt", "lin20.npt",
                   "linInf.1pt", "linInf.npt")
nb <- c(ode.1pt = grow(objs$ode, one_pt),
        ode.npt = grow(objs$ode, many_pt),
        lin20.1pt = grow(objs$lin20, one_pt),
        lin20.npt = grow(objs$lin20, many_pt),
        linInf.1pt = grow(objs$linInf, one_pt),
        linInf.npt = grow(objs$linInf, many_pt))
for (r in seq_len(ROUNDS)) {
  res[r, "ode.1pt"] <- one_pt(objs$ode, nb[["ode.1pt"]]) /
    nb[["ode.1pt"]]
  res[r, "ode.npt"] <- many_pt(objs$ode, nb[["ode.npt"]]) /
    nb[["ode.npt"]]
  res[r, "lin20.1pt"] <- one_pt(objs$lin20, nb[["lin20.1pt"]]) /
    nb[["lin20.1pt"]]
  res[r, "lin20.npt"] <- many_pt(objs$lin20, nb[["lin20.npt"]]) /
    nb[["lin20.npt"]]
  res[r, "linInf.1pt"] <- one_pt(objs$linInf, nb[["linInf.1pt"]]) /
    nb[["linInf.1pt"]]
  res[r, "linInf.npt"] <- many_pt(objs$linInf, nb[["linInf.npt"]]) /
    nb[["linInf.npt"]]
}
best <- apply(res, 2, min)
cat("block sizes:", paste(names(nb), nb, collapse = "  "), "\n")
cat("\nseconds per gradient, minimum over", ROUNDS, "rounds:\n")
print(signif(best, 4))
cat("\nrepeating ONE point understates the gradient by:\n")
for (a in c("ode", "lin20", "linInf")) {
  cat(sprintf("  %-8s %.2fx\n", a,
              best[[paste0(a, ".npt")]] / best[[paste0(a, ".1pt")]]))
}
cat("round spread (max/min):\n")
print(signif(apply(res, 2, max) / apply(res, 2, min), 4))

cat("\n--- where the whole fit's time goes ---\n")
for (a in names(forms)) {
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(frm(forms[[a]] + gaussian(), data = d,
                              start = ST, se = TRUE))
  tot <- proc.time()[["elapsed"]] - t0
  ev <- fit$opt$evaluations
  cat(sprintf("%-8s fit %8.2f s  iterations %3d  fn %4d  gr %4d",
              a, tot, fit$opt$iterations, ev[["function"]],
              ev[["gradient"]]))
  cat(sprintf("  logLik %.9f  code %d\n",
              as.numeric(stats::logLik(fit)), fit$opt$convergence))
}
