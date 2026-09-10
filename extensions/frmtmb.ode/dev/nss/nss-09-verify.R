# The change, through frm_ode() itself rather than a stand-in.
#
# A: the compartment grid against the exact limit, which frm_lincmt()
#    at its default n_ss = Inf supplies in closed form.
# B: the constructions that made each guard necessary, run through
#    frm_ode() rather than through the sweep's stand-in.
# C: ss_extrapolate = FALSE against the reference build of the base
#    commit, which must agree bit for bit. Run this script with
#    NSS_ARM=ref to write the reference side.
# D: the solve count, which is the load-independent cost.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-09-verify.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss"
ARM <- Sys.getenv("NSS_ARM", "lane")

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}
grid <- list(
  c(ke = 0.2,  k12 = 0.4,  k21 = 0.1,   ka = 1.1, ii = 8),
  c(ke = 0.15, k12 = 0.3,  k21 = 0.02,  ka = 1.0, ii = 24),
  c(ke = 0.15, k12 = 0.3,  k21 = 0.02,  ka = 1.0, ii = 12),
  c(ke = 0.1,  k12 = 0.2,  k21 = 0.008, ka = 1.0, ii = 24),
  c(ke = 0.08, k12 = 0.15, k21 = 0.003, ka = 1.0, ii = 24))
lamz <- function(g) {
  b <- g[["ke"]] + g[["k12"]] + g[["k21"]]
  (b - sqrt(b * b - 4 * g[["ke"]] * g[["k21"]])) / 2
}

ode_conc <- function(g, tt, ext, n = 20L) {
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = g[["ii"]],
                   ss = TRUE)
  z <- frm_ode(two_oral, init = list(0, 0, 0), times = tt,
               parms = list(g[["ke"]], g[["k12"]], g[["k21"]],
                            g[["ka"]]),
               events = ev, output = 2L, n_ss = n, ss_tol = 1,
               ss_extrapolate = ext)
  as.numeric(z) / 10
}
lin_conc <- function(g, tt) {
  ev <- data.frame(time = 0, state = "depot", value = 100,
                   ii = g[["ii"]], addl = 0L, ss = TRUE)
  as.numeric(frm_lincmt(parms = list(ke = g[["ke"]], k12 = g[["k12"]],
                                     k21 = g[["k21"]], ka = g[["ka"]],
                                     V = 10),
                        times = tt, ncmt = 2, depot = TRUE,
                        events = ev))
}

cat("\n=== A. the compartment grid, worst over one dosing interval ===\n")
cat(sprintf("%9s %5s %12s %12s %10s\n", "t_half_z", "ii", "truncated",
            "extrapolated", "gain"))
tabA <- NULL
for (g in grid) {
  tt <- seq(0, g[["ii"]], length.out = 25)
  ref <- lin_conc(g, tt)
  sc <- max(abs(ref))
  a <- max(abs(ode_conc(g, tt, FALSE) - ref)) / sc
  b <- max(abs(ode_conc(g, tt, TRUE) - ref)) / sc
  cat(sprintf("%9.1f %5g %12.3e %12.3e %10.2e\n",
              log(2) / lamz(g), g[["ii"]], a, b, a / b))
  tabA <- rbind(tabA, data.frame(t_half = log(2) / lamz(g),
                                 ii = g[["ii"]], trunc = a, ext = b))
}

cat("\n=== B. the constructions each guard exists for ===\n")
run_case <- function(tag, dyn, ini, parms, ii, tt, out, nref = 400L) {
  ev <- function(...) data.frame(time = 0, state = 1L, value = 100,
                                 ii = ii, ss = TRUE)
  f <- function(ext, n) {
    as.numeric(frm_ode(dyn, init = ini, times = tt, parms = parms,
                       events = ev(), output = out, n_ss = n,
                       ss_tol = 1, ss_extrapolate = ext))
  }
  ref <- f(FALSE, nref)
  a <- max(abs(f(FALSE, 20L) - ref)) / max(abs(ref))
  b <- max(abs(f(TRUE, 20L) - ref)) / max(abs(ref))
  cat(sprintf("  %-34s truncated %10.3e  extrapolated %10.3e\n",
              tag, a, b))
  invisible(c(a, b))
}
auc <- function(t, y, p) list(c(-p[1L] * y[1L],
                                p[1L] * y[1L] - p[2L] * y[2L],
                                y[2L] / 10))
one <- function(t, y, p) list(c(-p[1L] * y[1L],
                                p[1L] * y[1L] - p[2L] * y[2L]))
mm <- function(t, y, p) {
  cc <- y[2L] / 10
  list(c(-p[1L] * y[1L], p[1L] * y[1L] - p[2L] * cc / (p[3L] + cc)))
}
osc <- function(t, y, p) list(c(y[2L], -p[1L] * y[1L] - p[2L] * y[2L]))
run_case("AUC state, reading central", auc, list(0, 0, 0),
         list(1, 0.1), 24, seq(0, 24, length.out = 13), 2L)
run_case("AUC state, reading the AUC", auc, list(0, 0, 0),
         list(1, 0.1), 24, seq(0, 24, length.out = 13), 3L)
run_case("converged, ke*ii = 48", one, list(0, 0), list(3, 2), 24,
         seq(0, 24, length.out = 13), 2L)
run_case("Michaelis-Menten, Km = 50", mm, list(0, 0), list(1, 8, 50),
         24, seq(0, 24, length.out = 13), 2L, 4000L)
run_case("oscillator, ii = 1", osc, list(0, 0), list(1, 0.1), 1,
         seq(0, 1, length.out = 13), 1L, 4000L)

cat("\n=== C. ss_extrapolate = FALSE against the base commit ===\n")
bc <- list()
for (g in grid) {
  tt <- seq(0, g[["ii"]], length.out = 25)
  key <- paste0("g", g[["ke"]], "_", g[["k21"]], "_", g[["ii"]])
  bc[[key]] <- ode_conc(g, tt, FALSE)
}
saveRDS(bc, file.path(OUT, paste0("nss-09-backcompat-", ARM, ".rds")))
cat("  wrote", paste0("nss-09-backcompat-", ARM, ".rds"), "with",
    length(bc), "vectors\n")
ref_file <- file.path(OUT, "nss-09-backcompat-ref.rds")
if (ARM != "ref" && file.exists(ref_file)) {
  old <- readRDS(ref_file)
  cat("  identical to the base commit:",
      identical(bc, old), "\n")
  for (k in names(bc)) {
    d <- max(abs(bc[[k]] - old[[k]]))
    cat(sprintf("    %-28s max abs diff %s  bitwise %s\n", k,
                format(d), identical(bc[[k]], old[[k]])))
  }
}

cat("\n=== D. the solve count, which no clock is needed for ===\n")
count_solves <- function(ext, n) {
  cnt <- 0L
  tr <- trace(deSolve::lsoda, print = FALSE,
              tracer = function() cnt <<- cnt + 1L)
  on.exit(untrace(deSolve::lsoda), add = TRUE)
  g <- grid[[2L]]
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = g[["ii"]],
                   ss = TRUE)
  frm_ode(two_oral, init = list(0, 0, 0), times = c(0, 12),
          parms = list(g[["ke"]], g[["k12"]], g[["k21"]], g[["ka"]]),
          events = ev, output = 2L, n_ss = n, ss_tol = 1,
          ss_extrapolate = ext)
  cnt
}
for (ext in c(FALSE, TRUE)) {
  for (n in c(5L, 20L)) {
    cat(sprintf("  ss_extrapolate = %-5s n_ss = %-3d solves %d\n",
                ext, n, count_solves(ext, n)))
  }
}
