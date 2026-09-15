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
r <- tryCatch(posterior_summary(f), error = function(e) e)
if (inherits(r, "error")) say("posterior_summary(frmtmb_fit) ERROR: %s\n",
                              conditionMessage(r))
else say("posterior_summary(frmtmb_fit) rows=%d cols=%s\n", nrow(r),
         paste(colnames(r), collapse = ","))
for (g in c("bayes_R2", "hypothesis", "pp_check", "variables", "ngrps",
            "prior_summary", "as_draws", "nvariables", "loo_compare")) {
  v <- tryCatch({ x <- do.call(g, list(f)); paste("ok", class(x)[1]) },
                error = function(e) paste("ERR", substr(conditionMessage(e), 1, 70)))
  say("  %-16s %s\n", g, v)
}
cat("GENREVDONE\n")
