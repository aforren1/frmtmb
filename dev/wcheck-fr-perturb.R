# Lane wcheck, fix round: the new guard-absent fits of the tests (the
# reviewer's R1 under a centered Intercept prior, the bounded censored
# mean at E 2 and 6) and the R2 refusal, on data scaled by
# 1 + 2e-15 * k, k = 0..24.
#   Rscript dev/wcheck-fr-perturb.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
fit <- function(...) {
  w <- character(); m <- character()
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  }), error = function(e) e)
  olre <- "gives every observation its own random effect"
  if (inherits(f, "error")) {
    return(if (grepl("has no maximum|Laplace approximation",
                     conditionMessage(f))) "REFUSED" else "ERROR")
  }
  sprintf("code %d, %d other warnings, %d messages", f$opt$convergence,
          sum(!grepl(olre, w, fixed = TRUE)), length(m))
}
P <- function(...) set_prior(..., dpar = "sigma")
set.seed(20261111)
dc0 <- data.frame(obs = factor(seq_len(80)), x = stats::rnorm(80),
                  fac = factor(rep(c("a", "b"), 40)))
dc0$y <- 0.5 + 0.4 * dc0$x +
  stats::rnorm(80, 0, ifelse(dc0$fac == "a", 0.6, 1.2))
dc0$w <- ifelse(dc0$fac == "a", 0.8, 1.05)
set.seed(20261111)
dz0 <- data.frame(obs = factor(seq_len(80)), x = stats::rnorm(80))
dz0$y <- 1 + 0.5 * dz0$x + stats::rnorm(80)
q8 <- unname(stats::quantile(dz0$y, 0.8))
dz0$cc <- ifelse(dz0$y > q8, "right", "none")
nu0 <- sum(dz0$cc == "none")
fc <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
fz <- bf(yc | weights(w) + cens(cc) ~ 1 + (1 | obs))
pz <- set_prior("uniform(-5, 0.5)", class = "Intercept", lb = -5, ub = 0.5)
out <- list()
for (k in 0:24) {
  s <- 1 + 2e-15 * k
  dc <- dc0; dc$y <- dc$y * s; dc$x <- dc$x * s
  dz <- dz0; dz$y <- dz$y * s; dz$x <- dz$x * s
  dz$yc <- pmin(dz$y, q8 * s)
  out[[k + 1L]] <- c(
    center_R1 = fit(fc, family = gaussian(), data = dc,
                    prior = P("normal(0, 1)", class = "Intercept")),
    center_R2 = fit(fc, family = gaussian(), data = dc),
    cens_unif_E2 = fit(fz, family = gaussian(), prior = pz,
                       data = transform(dz, w = ifelse(cc == "none",
                                                       1 + 2 / nu0, 1))),
    cens_unif_E6 = fit(fz, family = gaussian(), prior = pz,
                       data = transform(dz, w = ifelse(cc == "none",
                                                       1 + 6 / nu0, 1))))
}
tab <- do.call(rbind, out)
for (j in colnames(tab)) {
  cat(sprintf("%-14s %s\n", j, paste(names(table(tab[, j])),
                                      table(tab[, j]), sep = " x",
                                      collapse = "; ")))
}
cat("DONE\n")
