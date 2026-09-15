# Lane eamhier: arm B, the two parameterizations PAIRED by seed.
#
# Arm B is the falsification design: `ndt` is the same for every subject
# and only the boundary varies, so the per-subject floors vary and the
# truth does not. dev/ndt-findings.md recorded, on 8 seeds at 15 x 250,
# a phantom spread of 6.40 ms under the per-group bound against 0.00
# under the global one, and a population `ndt` 23.58 ms low under the
# per-group bound against 23.87 ms low under the global one. Item 2.1
# asks whether that matters at the designs the field produces.
#
# The two arms are fitted on the SAME data, so every comparison here is
# paired and the paired standard error is the one to read.
#
# Run:  Rscript --vanilla dev/eamhier-scripts/eamhier-summary-B.R <dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
source("dev/eamhier-scripts/eamhier-summary-read.R")

d <- read_records(args)
d <- d[d$status == "ok", , drop = FALSE]
pg <- d[d$bound == "pg", , drop = FALSE]
gl <- d[d$bound == "gl", , drop = FALSE]
j <- match(pg$seed, gl$seed)
ok <- !is.na(j)
pg <- pg[ok, , drop = FALSE]
gl <- gl[j[ok], , drop = FALSE]
cat("arm B, ", pg$ns[1L], " subjects x ", pg$nt[1L], " trials, ",
    nrow(pg), " paired seeds (", min(pg$seed), " to ", max(pg$seed),
    ")\n", sep = "")

pair <- function(tag, a, b, unit = "") {
  dd <- a - b
  n <- length(dd)
  cat(sprintf(
    paste0("  %-24s pg %9.4f (sd %7.4f)  gl %9.4f (sd %7.4f)",
           "  diff %8.4f +- %6.4f %s\n"),
    tag, mean(a), stats::sd(a), mean(b), stats::sd(b), mean(dd),
    stats::sd(dd) / sqrt(n), unit))
}

pair("population ndt bias, ms", 1000 * (pg$ndt_hat - 0.25),
     1000 * (gl$ndt_hat - 0.25))
pair("reported sd(ndt), ms", pg$ndt_sd_hat_ms, gl$ndt_sd_hat_ms)
pair("per-subject rmse, ms", pg$ndt_rmse_ms, gl$ndt_rmse_ms)
pair("logLik", pg$logLik, gl$logLik)
pair("sd(floor), ms", pg$floor_sd_ms, gl$floor_sd_ms)
cat(sprintf("  per-group better on %d of %d seeds\n",
            sum(pg$logLik > gl$logLik), nrow(pg)))
cat(sprintf(
  "  ndt component under 1e-4: pg %d, gl %d of %d\n",
  sum(pg$sd2 < 1e-4, na.rm = TRUE),
  sum(gl$sd2 < 1e-4, na.rm = TRUE), nrow(pg)))
cat(sprintf("  convergence 0: pg %d, gl %d | pdHess TRUE: pg %d, gl %d\n",
            sum(pg$conv == 0), sum(gl$conv == 0),
            sum(pg$pdHess == "TRUE"), sum(gl$pdHess == "TRUE")))

# The phantom spread has a prediction attached to it: it is the FLOORS
# and nothing else, so it should track sd(floor) times the fitted
# fraction. That is checkable per replicate rather than asserted.
cat(sprintf("  correlation of pg sd(ndt) with sd(floor): %.4f\n",
            stats::cor(pg$ndt_sd_hat_ms, pg$floor_sd_ms)))
cat(sprintf("  pg sd(ndt) / (fraction x sd(floor)): mean %.4f, sd %.4f\n",
            mean(pg$ndt_sd_hat_ms / (pg$ndt_frac * pg$floor_sd_ms)),
            stats::sd(pg$ndt_sd_hat_ms /
                        (pg$ndt_frac * pg$floor_sd_ms))))
