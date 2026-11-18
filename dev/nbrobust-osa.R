# Lane nbrobust: the routes of nb_lpmf() other than a data count: the
# one-step residuals with the count on the tape ("oneStepGeneric"), the
# default "cdf" route and simulate(), on base and lane.
#   Rscript dev/nbrobust-osa.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(3)
d <- data.frame(x = stats::rnorm(60))
d$y <- stats::rnbinom(60, size = 2, mu = exp(1 + 0.5 * d$x))
d$y[c(3, 7)] <- c(31L, 45L)
fams <- list(negbinomial = negbinomial(), nbinom1 = nbinom1(),
             geometric = geometric(),
             zinb = zero_inflated_negbinomial(),
             hnb = hurdle_negbinomial(),
             nb_identity = negbinomial(link = "identity"))
for (nm in names(fams)) {
  fit <- suppressWarnings(suppressMessages(
    frm(y ~ x, family = fams[[nm]], data = d)))
  cat(sprintf("%s logLik %.12f code %d\n", nm, as.numeric(logLik(fit)),
              fit$opt$convergence))
  for (m in c("cdf", "oneStepGeneric")) {
    a <- list(fit, type = "osa", osa_method = m)
    if (m != "cdf") a$discreteSupport <- 0:120
    r <- tryCatch(suppressWarnings(do.call(residuals, a)),
                  error = function(e) paste("ERROR", conditionMessage(e)))
    if (is.character(r)) cat("  osa", m, r, "\n") else
      cat(sprintf("  osa %s sum %.10f sumsq %.10f NA %d\n", m,
                  sum(r, na.rm = TRUE), sum(r^2, na.rm = TRUE),
                  sum(is.na(r))))
  }
  s <- simulate(fit, nsim = 2, seed = 1)
  cat("  simulate head", paste(unlist(s)[1:6], collapse = " "), "\n")
}
cat("END\n")
