# Reviewer, lane wt-conditions: code in the installed interop packages
# that names simpleError, simpleWarning, simpleMessage or
# simpleCondition, which a classed frmtmb condition no longer matches.
#   Rscript dev/conditions-rev-interop-grep.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkgs <- c("emmeans", "insight", "marginaleffects", "broom.mixed", "broom",
          "testthat", "parameters", "performance", "bayestestR",
          "modelbased", "effectsize", "datawizard", "ggeffects", "lme4",
          "glmmTMB", "loo", "posterior", "bridgesampling", "brms",
          "car", "multcomp", "sandwich", "lmtest", "MuMIn", "DHARMa",
          "see", "report", "easystats", "knitr", "rmarkdown", "evaluate",
          "rlang", "purrr", "future", "furrr", "mice", "tidybayes",
          "shinystan", "bayesplot", "projpred", "rstan", "cmdstanr")
pat <- "simpleError|simpleWarning|simpleMessage|simpleCondition"
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("%-16s not installed\n", p)); next
  }
  ns <- asNamespace(p)
  hits <- character()
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns)
    if (!is.function(f) || is.primitive(f)) next
    txt <- deparse(f)
    i <- grep(pat, txt)
    if (length(i)) hits <- c(hits, paste0(nm, ": ", trimws(txt[i])))
  }
  cat(sprintf("%-16s %s  hits %d\n", p, format(packageVersion(p)),
              length(hits)))
  if (length(hits)) cat(paste0("    ", hits), sep = "\n")
}
