# Review of lane priors2, item 3: class theta on a model without a
# mixture, through frm() and frm_sample(), not only validate_prior().
# Rscript priors2-review-theta2.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(5)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, 15)),
                tf = factor(rep(1:3, 50)))
d$y2 <- rnorm(n)
try_msg <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:",
                                                conditionMessage(e)))
  if (!is.character(r)) {
    r <- paste("accepted; code", r$opt$convergence %||% "NA",
               "n entries", length(r$prior_entries %||% list()))
  }
  cat(lab, ":", substr(r, 1, 300), "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
ft <- bf(y2 ~ x + (1 | g))
try_msg("frm full, theta empty, no mixture",
        frm(ft, data = d, prior = set_prior("", class = "theta")))
try_msg("frm full, theta dirichlet(2), no mixture",
        frm(ft, data = d, prior = set_prior("dirichlet(2)", class = "theta")))
try_msg("frm dry_run objective, theta dirichlet(2), no mixture",
        frm(ft, data = d, prior = set_prior("dirichlet(2)", class = "theta"),
            dry_run = "objective"))
try_msg("frm dry_run objective, normal b (control)",
        frm(ft, data = d, prior = set_prior("normal(0, 1)", class = "b",
                                            coef = "zz"),
            dry_run = "objective"))
try_msg("frm full, empty b coef zz (control)",
        frm(ft, data = d, prior = set_prior("", class = "b", coef = "zz")))
try_msg("frm full, empty simo (control, no mo)",
        frm(ft, data = d, prior = set_prior("", class = "simo", coef = "mox1")))
try_msg("frm full, empty thetacov, no toep (control)",
        frm(ft, data = d, prior = set_prior("", class = "thetacov")))
