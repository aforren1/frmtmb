LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...))
set.seed(2026)
n <- 400L
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))
f <- frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd)
say("LIB %s\n", LIB)
for (g in c("bayes_R2", "pp_check", "variables", "ngrps",
            "prior_summary", "as_draws", "nvariables", "loo_compare",
            "expose_functions", "LOO", "WAIC", "loo", "waic")) {
  v <- tryCatch({ x <- do.call(g, list(f))
                  paste("ok", class(x)[1],
                        if (is.numeric(x)) sprintf("[%.4g..%.4g]", min(x), max(x)) else "") },
                error = function(e) paste("ERR", substr(conditionMessage(e), 1, 60)))
  say("  %-18s %s\n", g, v)
}
h <- tryCatch(hypothesis(f, "x = 0"), error = function(e) e)
if (inherits(h, "error")) {
  say("  hypothesis         ERR %s\n", substr(conditionMessage(h), 1, 60))
} else {
  say("  hypothesis         Estimate %.6f (link scale would be ~0.4)\n",
      h$hypothesis$Estimate[1])
}
say("  bayes_R2 value: %s\n",
    paste(signif(tryCatch(bayes_R2(f), error = function(e) NA), 6), collapse = ","))
cat("GENREVDONE\n")
