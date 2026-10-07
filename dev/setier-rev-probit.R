# Reviewer of lane setier: complete separation with a probit link is not
# named (dev/setier-rev-sep2.R). Why: the fitted tails, the certificate.
#   Rscript dev/setier-rev-probit.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
set.seed(4)
d <- data.frame(x = rnorm(60), z = rnorm(60))
d$yi <- as.integer(d$x > 0)
for (lk in c("probit", "logit")) {
  w <- character()
  f <- withCallingHandlers(frm(yi ~ x + z, family = bernoulli(lk), data = d),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  lp <- f$frame$linpreds[[1]]
  eta <- drop(as.matrix(lp$X) %*% f$estimates$beta[lp$idx])
  mu <- lp$link$linkinv(eta)
  tail <- ifelse(d$yi == 1, 1 - mu, mu)
  cat(lk, ": code", f$opt$convergence, f$opt$message, "| beta",
      signif(f$estimates$beta, 4), "\n  min |eta|", signif(min(abs(eta)), 3),
      "| tails: max", signif(max(tail), 3), " n <= 1e-3:", sum(tail <= 1e-3),
      " n <= 1e-9:", sum(tail <= 1e-9), " n == 0:", sum(tail == 0), "\n")
  cat("  link name:", lp$link$name, "| in sep_links:",
      lp$link$name %in% ns$sep_links, "\n")
  cat("  separation_check:", !is.null(ns$separation_check(f)), "\n")
  cat("  W:", substr(w, 1, 150), "\n")
}
