# The Phase 0 ode design, before and after, in ONE process.
#
# The design, the seed and the starting values are
# extensions/frmtmb.ode/tests/testthat/test-scale.R's, unchanged: 100
# subjects x 8 samples, a depot and a central state, twice-daily dosing
# for seven days written as one `ss` row plus `ii`/`addl`, n_ss = 20,
# `frm(se = TRUE)`. `dev/scale-findings.md` records 3948 s for the
# frm_ode() arm.
#
# Usage: Rscript lincmt-phase0.R <arms> <n_subject>
#   arms: "ode", "lin", or "both"
args <- commandArgs(trailingOnly = TRUE)
ARMS <- if (length(args) >= 1L) args[[1]] else "both"
NS <- if (length(args) >= 2L) as.integer(args[[2]]) else 100L

.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB); library(frmtmb); library(frmtmb.ode)
})
cat("frmtmb.ode", as.character(packageVersion("frmtmb.ode")), "from",
    dirname(getNamespaceInfo("frmtmb.ode", "path")), "\n")
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")
cat("started", format(Sys.time()), "\n")

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
f_lin <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses, n_ss = 20L),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
f_lin_inf <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}

run <- function(label, form) {
  cat("\n=====", label, "=====\n")
  t0 <- proc.time()[["elapsed"]]
  fit <- frm(form + gaussian(), data = d, start = ST, se = TRUE)
  el <- proc.time()[["elapsed"]] - t0
  b <- unlist(fixef(fit))
  ci <- suppressWarnings(stats::confint(fit))
  j <- grep("lke", rownames(ci), fixed = TRUE)
  vc <- VarCorr(fit)
  cat("wall clock (s):", format(el, digits = 6), "\n")
  cat("logLik:", format(as.numeric(stats::logLik(fit)), digits = 12),
      "\n")
  cat("ka:", exp(unname(b["lka.(Intercept)"])), " true", TR$ka, "\n")
  cat("ke:", exp(unname(b["lke.(Intercept)"])), " true", TR$ke,
      " CI", exp(as.numeric(ci[j[1L], 1:2])), "\n")
  cat("V :", exp(unname(b["lV.(Intercept)"])), " true", TR$V, "\n")
  cat("sd_lka:", sqrt(vc[[1L]][1L, 1L]), " sd_lke:",
      sqrt(vc[[2L]][1L, 1L]), "\n")
  cat("nlminb code:", fit$opt$convergence,
      " max |grad|:", format(max(abs(fit$obj$gr(fit$opt$par)))), "\n")
  cat("finished", format(Sys.time()), "\n")
  invisible(list(el = el, ll = as.numeric(stats::logLik(fit)),
                 par = fit$opt$par))
}

out <- list()
if (ARMS %in% c("lin", "both")) {
  out$lin <- run("AFTER: frm_lincmt(), n_ss = 20", f_lin())
  out$inf <- run("AFTER: frm_lincmt(), n_ss = Inf", f_lin_inf())
}
if (ARMS %in% c("ode", "both")) {
  out$ode <- run("BEFORE: frm_ode()", f_ode())
}
if (!is.null(out$ode) && !is.null(out$lin)) {
  cat("\n===== summary =====\n")
  cat("before (frm_ode)      ", format(out$ode$el, digits = 6), "s\n")
  cat("after  (n_ss = 20)    ", format(out$lin$el, digits = 6), "s\n")
  cat("after  (n_ss = Inf)   ", format(out$inf$el, digits = 6), "s\n")
  cat("speedup n_ss = 20     ",
      format(out$ode$el / out$lin$el, digits = 5), "x\n")
  cat("speedup n_ss = Inf    ",
      format(out$ode$el / out$inf$el, digits = 5), "x\n")
  cat("logLik difference     ",
      format(out$ode$ll - out$lin$ll, digits = 6), "\n")
  cat("max |par difference|  ",
      format(max(abs(out$ode$par - out$lin$par))), "\n")
}
cat("\ntotal session", format(proc.time()[["elapsed"]], digits = 6),
    "s\n")
