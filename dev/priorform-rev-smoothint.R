# Reviewer, lane wt-priorform: brms refuses x:s(z) and x:gp(z) as
# invalid terms. What do base and lane frmtmb build for them?
#   Rscript dev/priorform-rev-smoothint.R brms|ref|lane      seed 20260916
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(switch(mode, ref = "C:/Users/adf44/source/r/rellib-r3",
                   lane = "C:/Users/adf44/source/r/priorform-lib", NULL),
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (mode == "brms") suppressMessages(library(brms)) else suppressMessages(library(frmtmb))
set.seed(20260916)
n <- 200
d <- data.frame(x = rnorm(n), z = runif(n, 0, 3))
d$y <- 1 + d$x * sin(d$z) + rnorm(n, 0, .3)
fs <- c("y ~ x:s(z)", "y ~ x * s(z)", "y ~ x:gp(z)", "y ~ x * gp(z)",
        "y ~ s(z, by = x)", "y ~ x + s(z)", "y ~ x:mo(z)", "y ~ I(s(z))")
for (f in fs) {
  r <- tryCatch({
    if (mode == "brms") {
      suppressMessages(stancode(as.formula(f), data = d))
      "accept"
    } else {
      fit <- suppressWarnings(suppressMessages(frm(bf(as.formula(f)) + gaussian(), data = d)))
      lp <- fit$frame$linpreds[[1]]
      sprintf("accept: X cols [%s]; smooths %d; gp %d; logLik %.3f",
              paste(colnames(lp$X), collapse = ", "), length(lp$smooth %||% list()),
              length(lp$gp %||% list()), as.numeric(logLik(fit)))
    }
  }, error = function(e) paste("refuse:", substr(gsub("\n", " ", conditionMessage(e)), 1, 120)))
  cat(sprintf("%-18s %s\n", f, r))
}
