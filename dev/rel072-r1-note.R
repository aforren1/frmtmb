# Consolidation item 1 (lane prilabel's review, R1): the note of
# frm_sample() for a specification that applies to no parameter, when a
# coef row of the user's own shadows their class row. Data as in
# test-default-priors-brms.R (seed 51).
#   Rscript dev/rel072-r1-note.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "from",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(51)
n <- 80
d <- data.frame(x = runif(n), g = factor(sample(1:8, n, TRUE)))
d$y <- sin(4 * d$x) + rnorm(8, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)
uf <- frm(bf(y ~ gp(x) + (1 | g)), d, dry_run = "objective")
invisible(withCallingHandlers(
  frmtmb.sample:::sample_resolve_priors(uf,
    set_prior("normal(0, 1)", class = "sd") +
      set_prior("normal(0, 2)", class = "sd", group = "g",
                coef = "Intercept")),
  message = function(m) {
    x <- strsplit(conditionMessage(m), "\n")[[1]]
    cat(grep("applies to no parameter", x, value = TRUE), sep = "\n")
    invokeRestart("muffleMessage")
  }))
