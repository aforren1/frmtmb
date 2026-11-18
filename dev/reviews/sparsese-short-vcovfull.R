# Short check of lane sparsese: vcov(full = TRUE) of a REML fit whose group
# sd is at its boundary (the indep design, 100 groups, group sd 0, seed 1).
# The boundary message says the sd has NaN in vcov(); is vcov(full = TRUE)'s
# theta row NaN or finite?
#   Rscript sparsese-short-vcovfull.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
ng <- 100L
set.seed(1)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, 0)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
for (REML in c(TRUE, FALSE)) {
  fit <- withCallingHandlers(frm(fo, data = d, REML = REML),
    warning = function(w) {
      cat("W:", substr(conditionMessage(w), 1, 120), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("M:", substr(conditionMessage(m), 1, 120), "\n")
      invokeRestart("muffleMessage")
    })
  cat("REML", REML, "theta", fit$opt$par[["theta"]], "lost:",
      names(frmtmb:::sdr_of(fit)$se_lost), "\n")
  print(suppressWarnings(vcov(fit, full = TRUE)))
  print(suppressWarnings(summary(fit))$random)
}
