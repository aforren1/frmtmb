# Fits of lane dblbar's `||` cases (data seed 20261013): estimates,
# layouts of VarCorr(), ranef(), coef() and summary(), predictions at
# new levels and under re_formula, and the same model spelled
# (x | gr(g, cor = FALSE)) and diag(x | g).
# Usage: Rscript dev/dblbar-fits.R base|lane
# Output: dev/dblbar-log/fits-<arm>.rds and fits-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
sink(file.path(wt, sprintf("dev/dblbar-log/fits-%s.txt", arm)), split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
nd <- data.frame(x = c(0.3, -0.5), z = c(0.2, 0.1), f = c("b", "c"),
                 g = c("NEW1", "NEW2"), g1 = "p1", hg = "u")
cases <- list(
  num1 = list(y ~ x + (1 + x || g), y ~ x + (1 + x | gr(g, cor = FALSE)),
              y ~ x + diag(1 + x | g), y ~ x + (1 | g) + (0 + x | g)),
  fac1 = list(y ~ f + (1 + f || g), y ~ f + (1 + f | gr(g, cor = FALSE)),
              y ~ f + diag(1 + f | g)),
  mix = list(y ~ z + f + (1 + z + f || g),
             y ~ z + f + (1 + z + f | gr(g, cor = FALSE)),
             y ~ z + f + diag(1 + z + f | g))
)
one <- function(form) {
  fit <- tryCatch(frm(bf(form), family = gaussian(), data = dd),
                  error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  pv <- frmtmb:::hyp_vals_only(fit)
  vals <- frmtmb:::hyp_env_vals(fit, pv$vals, pv$comp)
  pred <- function(...) {
    tryCatch(fitted(fit, newdata = nd, allow_new_levels = TRUE, ...),
             error = function(e) conditionMessage(e))
  }
  set.seed(3)
  p_unc <- pred(sample_new_levels = "uncertainty")
  set.seed(3)
  p_gau <- pred(sample_new_levels = "gaussian")
  p_pop <- tryCatch(fitted(fit, newdata = nd, re_formula = NA),
                    error = function(e) conditionMessage(e))
  p_int <- tryCatch(fitted(fit, newdata = dd[1:6, ], re_formula = ~ (1 | g)),
                    error = function(e) conditionMessage(e))
  p_all <- tryCatch(fitted(fit, newdata = dd[1:6, ],
                           re_formula = stats::as.formula(
                             paste("~", deparse1(form[[3]][[3]])))),
                    error = function(e) conditionMessage(e))
  list(logLik = as.numeric(logLik(fit)), opt = fit$opt$par,
       vals = unlist(vals), theta = fit$estimates$theta,
       vc = VarCorr(fit), ranef_dimnames = lapply(ranef(fit), dimnames),
       coef_dimnames = lapply(coef(fit), dimnames),
       summary_random = lapply(summary(fit)$random, rownames),
       blocks = vapply(fit$frame$re_blocks, function(b) {
         paste0(b$covstruct, "[", paste(b$cnms, collapse = ","), "]")
       }, ""),
       p_unc = p_unc, p_gau = p_gau, p_pop = p_pop, p_int = p_int,
       p_all = p_all, fitted_in = fitted(fit)[1:6, 1],
       ranef_g = tryCatch(ranef(fit)$g, error = function(e) conditionMessage(e)))
}
res <- list()
for (nm in names(cases)) {
  res[[nm]] <- lapply(cases[[nm]], one)
  names(res[[nm]]) <- vapply(cases[[nm]], deparse1, "")
  for (k in seq_along(res[[nm]])) {
    r <- res[[nm]][[k]]
    cat("\n==== ", names(res[[nm]])[k], "\n", sep = "")
    if (!is.null(r$error)) {
      cat("ERROR:", r$error, "\n")
      next
    }
    cat("logLik", format(r$logLik, digits = 15), " blocks",
        paste(r$blocks, collapse = " "), "\n")
    print(signif(r$vals, 8))
    cat("summary random:", paste(unlist(r$summary_random), collapse = " "),
        "\n")
    cat("VarCorr elements:", paste(names(r$vc$g), collapse = " "), "\n")
  }
}
sink()
saveRDS(res, file.path(wt, sprintf("dev/dblbar-log/fits-%s.rds", arm)))
