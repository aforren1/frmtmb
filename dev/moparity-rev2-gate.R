# Reviewer, re-check (m3): the new test-reparam.R correlated gate
# (|z| < 5 in MCSE units on every mean and sd, centered vs
# non-centered), its false-alarm rate over more seeds, and its power
# against a wiring bug: the non-centered factor transposed, b = L' z,
# injected in this process only (the registry entry of `us`).
#   Rscript dev/moparity-rev2-gate.R <design> <mutate 0|1> first last
# design: orig (12 groups, 120 rows, sds 0.7 0.4, independent),
#         new (30, 240, 0.8 0.6, independent),
#         newcor (30, 240, 0.8 0.6, correlation 0.6)
args <- commandArgs(TRUE)
design <- args[1]
mutate <- identical(args[2], "1")
seeds <- as.integer(args[3]):as.integer(args[4])
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("design", design, "mutate", mutate, "frmtmb", find.package("frmtmb"), "\n")
if (mutate) {
  ns <- asNamespace("frmtmb")
  reg <- get("covstruct_registry", ns)
  orig <- reg[["us"]]$chol_L
  reg[["us"]]$chol_L <- function(theta, blk) t(orig(theta, blk))
  unlockBinding("covstruct_registry", ns)
  assign("covstruct_registry", reg, envir = ns)
  lockBinding("covstruct_registry", ns)
}
set.seed(31)
des <- switch(design, orig = c(12, 120, 0.7, 0.4, 0), new = c(30, 240, 0.8, 0.6, 0),
              newcor = c(30, 240, 0.8, 0.6, 0.6))
ng <- des[1]; nn <- des[2]
dd <- data.frame(x = stats::rnorm(nn),
                 g = factor(rep(seq_len(ng), length.out = nn)))
u <- matrix(stats::rnorm(2 * ng), 2)
u[2, ] <- des[5] * u[1, ] + sqrt(1 - des[5]^2) * u[2, ]
dd$y <- stats::rnorm(nn, 1 + 0.5 * dd$x + des[3] * u[1, dd$g] +
                       des[4] * u[2, dd$g] * dd$x, 1)
q <- function(f, v) suppressWarnings(f(v))
for (s in seeds) {
  mk <- function(ncp) {
    suppressWarnings(suppressMessages(
      frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
                 chains = 1, iter = 1000, refresh = 0, seed = s,
                 reparameterize = ncp)))
  }
  cen <- mk(FALSE)
  ncp <- mk(TRUE)
  zs <- c()
  for (nm in setdiff(colnames(cen$draws), "lp__")) {
    a <- cen$draws[, nm]; b <- ncp$draws[, nm]
    zm <- (mean(a) - mean(b)) / sqrt(q(posterior::mcse_mean, a)^2 + q(posterior::mcse_mean, b)^2)
    zsd <- (sd(a) - sd(b)) / sqrt(q(posterior::mcse_sd, a)^2 + q(posterior::mcse_sd, b)^2)
    zs[paste(nm, "mean")] <- zm
    zs[paste(nm, "sd")] <- zsd
  }
  w <- which.max(abs(zs))
  cat(sprintf("RES design %s mutate %d seed %d | blocks ncp %s | max |z| %.2f (%s) | pass %s | min ESS cen %.0f\n",
              design, mutate, s, format(ncp$reparam$blocks %||% NA), abs(zs[w]), names(zs)[w],
              all(abs(zs) < 5),
              min(apply(cen$draws[, setdiff(colnames(cen$draws), "lp__")], 2,
                        function(v) q(posterior::ess_bulk, v)))))
}
