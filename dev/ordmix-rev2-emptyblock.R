# Reviewer of lane ordmix (re-check, copy of dev/ordmix-rev-emptyblock.R): the draws_to_natural() empty-block idiom the
# lane found by reading (dev/ordmix-findings.md, "Not done"). An ordinal
# block with no internal parameter: sum-to-zero thresholds on a binary
# response, one threshold held at 0. Does frm_sample() store brms's
# b_Intercept[1] column? And the inverse, driven by hand on a matrix
# that does carry the name. Usage: Rscript ... <lane|base>. Seed 11.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1L + (0.8 * d$x + rlogis(n) > 0)
fit <- frm(bf(y ~ x), family = cumulative(threshold = "sum_to_zero"),
           data = d)
cat("fixef rows:", paste(rownames(fixef(fit)), collapse = " "), "\n")
cat("template tau_raw length:", length(fit$frame$par_template$tau_raw), "\n")
nc <- frmtmb.sample:::draws_natural_cols(fit)
for (o in nc$ordinal) {
  cat("block internal:", length(o$internal), " names:",
      paste(o$names, collapse = " "), "\n")
}
ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 200, refresh = 0,
                                  seed = 3,
                                  prior = set_prior("normal(0, 2)",
                                                    class = "b")))
cat("draw columns:", paste(colnames(as_draws_matrix(ds)), collapse = " "),
    "\n")
# the inverse on a matrix that carries the names brms would store
m <- cbind(b_x = c(0.5, 0.6), "b_Intercept[1]" = c(0, 0), lp__ = c(-1, -2))
inv <- tryCatch(frmtmb.sample:::draws_to_natural(m, fit, inverse = TRUE),
                error = function(e) paste("ERROR", conditionMessage(e)))
cat("inverse on b_x, b_Intercept[1], lp__ gives columns:",
    if (is.matrix(inv)) paste(colnames(inv), collapse = " ") else inv,
    " ncol", if (is.matrix(inv)) ncol(inv) else NA, "\n")
