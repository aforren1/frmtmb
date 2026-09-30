# frmtmb.sample on disc and threshold-structure fits: draw names, the
# round trip of the draws to the internal scale, the disc default prior,
# and the draws methods. Seed 20260930 for the data, 3 for the sampler.
# Output: dev/ordinal-log-sample-smoke.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb.sample from", find.package("frmtmb.sample"), "frmtmb from",
    find.package("frmtmb"), "\n")
set.seed(20260930)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
run <- function(lab, fit) {
  cat("\n==", lab, "==\n")
  msg <- character()
  ds <- withCallingHandlers(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }, warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  cat(msg, sep = "")
  cat("variables:", paste(variables(ds), collapse = " "), "\n")
  # the stored draws map back to the internal scale and forward again
  m <- ds$draws
  back <- frmtmb.sample:::draws_to_natural(
    frmtmb.sample:::draws_internal_matrix(ds), ds$fit)
  cat("round trip max rel diff:",
      max(abs(back - m[, colnames(back)]) / pmax(1, abs(m[, colnames(back)]))),
      " same columns:", identical(colnames(back), colnames(m)), "\n")
  print(fixef(ds))
  print(dim(posterior_epred(ds, ndraws = 5)))
  print(dim(log_lik(ds, ndraws = 5)))
  print(table(posterior_predict(ds, ndraws = 5)))
  print(prior_summary(ds))
  invisible(ds)
}
ds1 <- run("cumulative equidistant, disc ~ z",
           frm(bf(y ~ x, disc ~ z), family = cumulative(threshold =
                                                        "equidistant"),
               data = d,
               prior = set_prior("normal(0, 1)", class = "Intercept",
                                 dpar = "disc")))
print(posterior_summary(ds1, variable = "delta"))
print(hypothesis(ds1, "delta > 0.5", class = NULL)$hypothesis)
ds2 <- run("sratio sum_to_zero", frm(y ~ x,
                                     family = sratio(threshold = "sum_to_zero"),
                                     data = d))
b <- as.matrix(ds2, variable = paste0("b_Intercept[", 1:4, "]"))
cat("sum-to-zero in the draws, max |row sum|:", max(abs(rowSums(b))), "\n")
ds3 <- run("hurdle flexible", frm(yh ~ x, family = hurdle_cumulative(),
                                  data = d))
