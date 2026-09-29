# REVIEW script 07, claim 4: are pp_check()'s loo_* types broken on
# EVERY model, including on the reference build?
#
#   Rscript dev/arcovsample-rev-07-ppcheck.R <lane|ref>
#
# Three models, so "every model" is a measurement and not a single
# case: a plain gaussian fit with nothing special, one with a random
# effect, and (lane only, because the reference refuses log_lik there) a
# cov = FALSE ARMA fit. Draws come from a fixed matrix, no sampler.

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
newx <- all(c("arma_cond_resp", "arma_cond_dpars") %in%
              getNamespaceExports("frmtmb"))
cat("ARM ", arm, " newexports=", newx, " bayesplot ",
    format(packageVersion("bayesplot")), " loo ",
    format(packageVersion("loo")), "\n", sep = "")
stopifnot(identical(newx, identical(arm, "lane")))

set.seed(31L)
dd <- data.frame(g = factor(rep(1:5, each = 6L)), t = rep(1:6, 5L))
n <- nrow(dd)
dd$x <- rnorm(n)
dd$y <- 0.6 + 0.5 * dd$x + rnorm(n, 0, 0.8)

mk <- function(fit, nd = 20L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(0, nd, length(lab), dimnames = list(NULL, lab))
  set.seed(77L)
  for (j in seq_along(lab)) {
    nm <- lab[j]
    m[, j] <- if (nm == "lp__") 0
    else if (grepl("^sigma", nm) || grepl("^sd_", nm)) {
      exp(rnorm(nd, log(0.8), 0.05))
    } else if (grepl("^thetaac", nm)) rnorm(nd, 0.4, 0.05)
    else rnorm(nd, 0.5, 0.1)
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

cases <- list(
  plain = frm(bf(y ~ x), family = gaussian(), data = dd,
              dry_run = "objective"),
  re = frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
           dry_run = "objective"),
  arma = frm(bf(y ~ x + ar(t, g)), family = gaussian(), data = dd,
             dry_run = "objective"))

types <- c("dens_overlay", "stat", "loo_pit_overlay", "loo_pit",
           "loo_intervals", "loo_ribbon")
for (cn in names(cases)) {
  ds <- mk(cases[[cn]])
  cat("\n---- model ", cn, "\n", sep = "")
  for (ty in types) {
    r <- tryCatch({
      p <- suppressWarnings(suppressMessages(pp_check(ds, type = ty,
                                                      ndraws = 10L)))
      paste0("OK ", paste(class(p), collapse = "/"))
    }, error = function(e) paste0(class(e)[1L], ": ",
                                  conditionMessage(e)))
    cat("  [", ty, "] ", r, "\n", sep = "")
  }
}
cat("\nDONE\n")
