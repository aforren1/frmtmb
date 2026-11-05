# The ordinary spelling of the scale defect: a calendar year as a covariate,
# uncentered. y ~ year + (1 | g), year = 2000 + 5 z. The same fit on the centered
# column is the reference (the model is the same model, and an affine map of the
# column leaves the slope, its SE, theta and the logLik unchanged).
# Usage: Rscript dev/review-impsamp-year.R [after|before]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- commandArgs(TRUE)
if (!length(arm)) arm <- "after"
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
options(width = 200)
rows <- list()
for (fam in c("poisson", "bernoulli")) {
  for (s in 1:6) {
    set.seed(s)
    G <- 15L
    g <- factor(rep(seq_len(G), each = 10))
    z <- rnorm(150)
    u <- rnorm(G, 0, 0.8)
    eta <- if (fam == "poisson") 0.5 + 0.3 * z + u[g] else 0.2 + 0.6 * z + u[g]
    y <- if (fam == "poisson") rpois(150, exp(eta)) else rbinom(150, 1, plogis(eta))
    family <- if (fam == "poisson") poisson() else bernoulli()
    d <- data.frame(y, year = round(2000 + 5 * z), g)
    d$yc <- d$year - mean(d$year)
    w <- character(0)
    t0 <- Sys.time()
    f <- withCallingHandlers(
      tryCatch(frm(bf(y ~ year + (1 | g)), family = family, data = d, importance = 300L),
               error = function(e) e),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    r <- suppressWarnings(tryCatch(frm(bf(y ~ yc + (1 | g)), family = family, data = d,
                                       importance = 300L), error = function(e) e))
    if (inherits(r, "error")) r <- suppressWarnings(frm(bf(y ~ yc + (1 | g)), family = family, data = d))
    if (inherits(f, "error")) {
      rows[[length(rows) + 1L]] <- data.frame(arm, fam, s, rounds = NA, maxgrad = NA,
        moved = NA, dslope_z = NA, dtheta_z = NA, se_slope = NA, se_theta = NA, dll = NA,
        mcse = NA, nw = NA, secs = secs, note = substr(conditionMessage(f), 1, 80))
      next
    }
    sf <- sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
    sr <- sqrt(diag(frmtmb:::sdr_of(r)$cov.fixed))
    rows[[length(rows) + 1L]] <- data.frame(
      arm, fam, s, rounds = f$importance$rounds, maxgrad = f$importance$grad,
      moved = f$importance$moved,
      dslope_z = (f$opt$par[[2]] - r$opt$par[[2]]) / sr[[2]],
      dtheta_z = (f$opt$par[[3]] - r$opt$par[[3]]) / sr[[3]],
      se_slope = sf[[2]] / sr[[2]], se_theta = sf[[3]] / sr[[3]],
      dll = as.numeric(logLik(f)) - as.numeric(logLik(r)), mcse = f$importance$mcse,
      nw = length(w), secs = secs, note = substr(paste(w, collapse = " | "), 1, 80))
  }
}
out <- do.call(rbind, rows)
print(format(out, digits = 3), row.names = FALSE)
write.table(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/impsamp-out/year-%s.tsv", arm),
  sep = "\t", quote = FALSE, row.names = FALSE)
