# Punch round 1: NEWS edits (m1, m6, B1, B2), done by exact replacement.
p <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/NEWS.md"
s <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(s, gregexpr(old, s, fixed = TRUE))) == 1L)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1(paste0("\n\n* **`mixture(order = \"mu\")` with ordinal components is brms's shared\n",
            "  thresholds** (below); with components that have a `mu` intercept it\n",
            "  stays refused."), "")
rep1(paste0("  `fixed_Intercept`, which a sum-to-zero component centers.\n"),
     paste0("  `fixed_Intercept`, which a sum-to-zero component centers (with\n",
            "  components that have a `mu` intercept, `order = \"mu\"` stays\n",
            "  refused).\n"))
rep1(paste0("  parameter vectors to at most 2.5 ulp on 20 mixture shapes\n",
            "  (`dev/ordmix-lpcheck.R`), and `fixef()` has brms's rows in brms's\n"),
     paste0("  parameter vectors to at most 2.5 ulp on the 20 mixture shapes of\n",
            "  `dev/ordmix-lpcheck.R`; a probit component whose latent distance\n",
            "  from a threshold passes about 38 has a `NaN` density where brms's\n",
            "  is finite (one point of the review's three-component probit, sratio\n",
            "  and acat mixture). `fixef()` has brms's rows in brms's\n"))
rep1(paste0("  likelihood has a flat direction by construction warns: no predictor\n",
            "  in any parameter, hurdle components whose `hu<k>` and mixing weights\n",
            "  have none, or shared thresholds with no `mu<k>` predictor.\n"),
     paste0("  likelihood has a flat direction by construction warns: no predictor\n",
            "  in any parameter, or hurdle components whose `hu<k>` and mixing\n",
            "  weights have none, unless a prior holds them; or shared thresholds\n",
            "  with no `mu<k>` predictor on components of one family, link and\n",
            "  `disc`. A fit whose component ends at a degenerate boundary (a step\n",
            "  function of its predictors, or two thresholds that are the same\n",
            "  number) warns and names it: maximum likelihood for an ordinal\n",
            "  mixture often has its supremum there (22 of 160 simulated\n",
            "  two-component fits; the warning fires on 20 of them and on none of\n",
            "  the other 138), so compare several starts, or hold the components\n",
            "  with priors.\n"))
rep1(paste0("## Bug fixes\n\n* **The draws of a model whose location"),
     paste0("## Bug fixes\n\n",
            "* **A gradient that is not finite at the reported optimum** passed\n",
            "  `check_convergence()` in silence, because its test compares a finite\n",
            "  number with `grad_tol`. nlminb reports X-convergence there (a\n",
            "  three-component ordinal mixture whose probit component saturated,\n",
            "  gradient `Inf`). The fit now warns and names the parameters.\n\n",
            "* **The draws of a model whose location"))
writeLines(s, p)
cat("done\n")
