# Review of lane priors2, item 2: the refusal names frm_sample() as the
# way out; does frm_sample() on the formula route take the same prior?
# Rscript priors2-review-sample.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priors2/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
dc <- data.frame(z = rnorm(120))
dc$y <- 1 + 0.3 * dc$z + rnorm(120, sd = 0.3)
r <- tryCatch(suppressMessages(frm_sample(
  y ~ z, data = dc,
  prior = set_prior("gamma(0.5, 1)", class = "b", lb = 0),
  chains = 1, iter = 300, seed = 1, refresh = 0)), error = function(e) e)
if (inherits(r, "error")) {
  cat("frm_sample: ERROR", conditionMessage(r), "\n")
} else {
  m <- as.matrix(r)
  cat("frm_sample: drew", nrow(m), "draws; b_z mean",
      format(mean(m[, "b_z"]), digits = 6), "\n")
}
# frm() on the same class-wide row: still refused?
r2 <- tryCatch(frm(y ~ z, data = dc, prior = set_prior("gamma(0.5, 1)",
                                                       class = "b", lb = 0)),
               error = function(e) e)
cat("frm class-wide:", if (inherits(r2, "error")) conditionMessage(r2) else
  paste("code", r2$opt$convergence), "\n")
