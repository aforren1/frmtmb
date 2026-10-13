# Reviewer of lane optrap, final check: responses where the MAD of the
# least-squares residuals is 0 or tiny next to their sd, fitted as
# gaussian with (1 | g) (20 groups of 10), seeds 1 to 20, against lmer:
#   zeros90  y = 10 * Bernoulli(0.1) + 0.01 x (90 percent of rows at
#            one value up to a tiny slope)
#   zeros80i y = 10 * Bernoulli(0.2), intercept only (MAD 0, the sd
#            fallback)
#   ties     y = round(N(0, 0.6)) + 0.3 x (many ties)
#   counts   y = Poisson(0.5) + 0.2 x
#   Rscript dev/optrap-rev3-mad.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
pu <- tryCatch(get("predictor_units", asNamespace("frmtmb")),
               error = function(e) NULL)
X <- list()
for (cs in c("zeros90", "zeros80i", "ties", "counts")) for (s in 1:20) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200))
  b <- rnorm(20, 0, 0.3)[d$g]
  fo <- y ~ x + (1 | g)
  if (cs == "zeros90") d$y <- 10 * rbinom(200, 1, 0.1) + 0.01 * d$x + b
  if (cs == "zeros80i") {
    d$y <- 10 * rbinom(200, 1, 0.2)
    fo <- y ~ 1 + (1 | g)
  }
  if (cs == "ties") d$y <- round(rnorm(200, 0, 0.6) + b) + 0.3 * d$x
  if (cs == "counts") d$y <- rpois(200, 0.5) + 0.2 * d$x
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  l4 <- suppressMessages(suppressWarnings(lmer(fo, data = d, REML = FALSE,
    control = lmerControl(check.conv.singular = "ignore"))))
  u <- if (is.null(pu)) "none" else {
    uu <- pu(f$frame)
    if (length(uu)) format(signif(uu[[1]], 3)) else "1"
  }
  X[[length(X) + 1L]] <- data.frame(
    case = cs, seed = s, code = f$opt$convergence,
    dll = as.numeric(logLik(f)) - as.numeric(logLik(l4)),
    nwarn = sum(grepl("did not report|maximum absolute gradient", w)),
    unit = u, sdy = sd(d$y))
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev3-log/mad-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
for (cs in unique(X$case)) {
  x <- X[X$case == cs, ]
  cat(sprintf("%s %-9s | units %s | sd(y) median %.2f | convergence-warned %d | code != 0 %d | short of lmer > 1e-4 %d (min %.3g)\n",
              lab, cs, paste(sort(unique(x$unit)), collapse = ","),
              median(x$sdy), sum(x$nwarn > 0), sum(x$code != 0),
              sum(x$dll < -1e-4), min(x$dll)))
}
