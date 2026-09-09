# Attack 5: the n_ss finding, and its severity.
#
# The lane measured ONE schedule and reported 3.96e-03. The severity
# question is not what that schedule does, it is what the shortfall is
# a function of, and whether the bias it leaves changes an estimate.
# Both are cheap to answer with the closed form, because n_ss = Inf is
# the exact limit and n_ss = 20 is exactly what frm_ode() simulates.
#
# Script path: dev/rev-lincmt-nss.R. Seeds recorded at each fit.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

cat("\n=== A. the shortfall as a function of lambda_z * ii ===\n")
cat("Two compartments, oral, one ss row. lambda_z is the slow",
    "disposition\neigenvalue; the shortfall is",
    "frm_lincmt(n_ss = 20) against n_ss = Inf,\nworst over one",
    "dosing interval, relative to the trajectory's scale.\n\n")
lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}
short <- function(ke, k12, k21, ka, ii, n) {
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = ii,
                   addl = 0L, ss = TRUE)
  ts <- seq(0, ii, length.out = 25)
  P <- list(ke = ke, k12 = k12, k21 = k21, ka = ka, V = 10)
  a <- frm_lincmt(parms = P, times = ts, ncmt = 2, depot = TRUE,
                  events = ev, n_ss = n)
  b <- frm_lincmt(parms = P, times = ts, ncmt = 2, depot = TRUE,
                  events = ev)
  max(abs(a - b)) / max(abs(b))
}
cat(sprintf("%10s %8s %10s %12s %12s %12s %12s\n", "t_half_z", "ii",
            "lam_z*ii", "n_ss=20", "n_ss=40", "n_ss=80", "n_ss=200"))
grid <- list(
  c(ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, ii = 8),
  c(ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, ii = 24),
  c(ke = 0.5, k12 = 1.0, k21 = 0.05, ka = 1.1, ii = 12),
  c(ke = 0.15, k12 = 0.3, k21 = 0.02, ka = 1.0, ii = 24),
  c(ke = 0.1, k12 = 0.2, k21 = 0.008, ka = 1.0, ii = 24),
  c(ke = 0.08, k12 = 0.15, k21 = 0.003, ka = 1.0, ii = 24),
  c(ke = 0.05, k12 = 0.1, k21 = 0.001, ka = 1.0, ii = 168))
for (g in grid) {
  lz <- lamz(g[["ke"]], g[["k12"]], g[["k21"]])
  cat(sprintf("%10.1f %8g %10.4f %12.3e %12.3e %12.3e %12.3e\n",
              log(2) / lz, g[["ii"]], lz * g[["ii"]],
              short(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]],
                    g[["ii"]], 20L),
              short(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]],
                    g[["ii"]], 40L),
              short(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]],
                    g[["ii"]], 80L),
              short(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]],
                    g[["ii"]], 200L)))
}

cat("\n=== B. one compartment, which is the easy case ===\n")
short1 <- function(ke, ka, ii, n) {
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = ii,
                   addl = 0L, ss = TRUE)
  ts <- seq(0, ii, length.out = 25)
  P <- list(ke = ke, ka = ka, V = 10)
  a <- frm_lincmt(parms = P, times = ts, ncmt = 1, depot = TRUE,
                  events = ev, n_ss = n)
  b <- frm_lincmt(parms = P, times = ts, ncmt = 1, depot = TRUE,
                  events = ev)
  max(abs(a - b)) / max(abs(b))
}
cat(sprintf("%10s %8s %10s %12s\n", "t_half", "ii", "ke*ii",
            "n_ss=20"))
for (kk in c(0.15, 0.05, 0.02, 0.01, 0.005)) {
  cat(sprintf("%10.1f %8g %10.4f %12.3e\n", log(2) / kk, 12,
              kk * 12, short1(kk, 1.0, 12, 20L)))
}

cat("\n=== C. does the bias move an estimate? ===\n")
cat("Simulate from the EXACT steady state, then fit twice: once with",
    "\nn_ss = Inf (the truth) and once with n_ss = 20 (what frm_ode()",
    "\nsimulates). Both arms are the closed form, so the only",
    "difference\nbetween them is the run-in truncation.\n\n")
fit_pair <- function(ke, k12, k21, ka, ii, V, ns, seed, sd_obs) {
  set.seed(seed)
  tt <- ii * c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1)
  d <- data.frame(id = factor(rep(seq_len(ns), each = length(tt))),
                  time = rep(tt, ns))
  i <- as.integer(d$id)
  lke <- log(ke) + rnorm(ns, 0, 0.2)
  lka <- log(ka) + rnorm(ns, 0, 0.3)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = ii,
                   addl = 0L, ss = TRUE)
  mu <- numeric(nrow(d))
  for (j in seq_len(ns)) {
    k <- which(i == j)
    mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = k12,
                                     k21 = k21, ka = exp(lka[[j]]),
                                     V = V),
                        times = d$time[k], ncmt = 2, depot = TRUE,
                        events = ev)
  }
  d$conc <- mu + rnorm(nrow(d), 0, sd_obs)
  doses <- ev
  mk <- function(nss) {
    bd <- if (is.finite(nss))
      bf(conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                        k21 = exp(lk21), ka = exp(lka),
                                        V = exp(lV)),
                           times = time, group = id, ncmt = 2,
                           depot = TRUE, events = doses, n_ss = 20L),
         lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
         lk21 ~ 1, lV ~ 1, nl = TRUE)
    else
      bf(conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                        k21 = exp(lk21), ka = exp(lka),
                                        V = exp(lV)),
                           times = time, group = id, ncmt = 2,
                           depot = TRUE, events = doses),
         lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
         lk21 ~ 1, lV ~ 1, nl = TRUE)
    frm(bd + gaussian(), data = d,
        start = list(beta = c(log(ke), log(ka), log(k12), log(k21),
                              log(V))))
  }
  list(inf = mk(Inf), n20 = mk(20L))
}
run <- function(tag, ke, k12, k21, ka, ii, V, seed) {
  z <- tryCatch(fit_pair(ke, k12, k21, ka, ii, V, 30L, seed, 0.15),
                error = function(e) conditionMessage(e))
  if (is.character(z)) { cat(tag, "ERROR:", substr(z, 1, 80), "\n")
    return(invisible(NULL)) }
  a <- fixef(z$inf); b <- fixef(z$n20)
  ea <- if (is.matrix(a)) a[, 1L] else unlist(a)
  eb <- if (is.matrix(b)) b[, 1L] else unlist(b)
  nm <- names(ea)
  if (is.null(nm)) nm <- paste0("beta", seq_along(ea))
  cat(tag, "  seed", seed, "\n")
  cat(sprintf("  %-16s %12s %12s %10s\n", "parameter", "n_ss = Inf",
              "n_ss = 20", "bias %"))
  for (j in seq_along(nm))
    cat(sprintf("  %-16s %12.5f %12.5f %10.2f\n", nm[[j]],
                ea[[j]], eb[[j]],
                100 * (eb[[j]] - ea[[j]]) / abs(ea[[j]])))
  cat(sprintf("  %-14s %12.4f %12.4f\n", "logLik",
              as.numeric(logLik(z$inf)), as.numeric(logLik(z$n20))))
}
run("slow peripheral, ii = 24", 0.1, 0.2, 0.008, 1.0, 24, 10, 101)
run("the package's own 2 cmt schedule, ii = 8", 0.2, 0.4, 0.1, 1.1, 8,
    10, 102)
