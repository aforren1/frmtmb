# Reviewer of lane optrap, item 8 widened: on y ~ x + z (gaussian, no
# random effect), where lm() gives the exact maximum, how far each build
# stops from it, in the fit's own standard errors, by the response's
# spread. Seeds 1 to 20.
#   Rscript dev/optrap-rev-lmprec.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
X <- list()
for (sdv in c(1, 1.5, 2, 4, 8, 100, 1e3, 1e6)) for (s in 1:20) {
  set.seed(s)
  n <- 100
  d <- data.frame(x = rnorm(n), z = runif(n, 0, 10))
  d$y <- (1 + 0.5 * d$x - 0.2 * d$z + rnorm(n)) * sdv
  f <- suppressWarnings(frm(bf(y ~ x + z), family = gaussian(), data = d))
  l <- lm(y ~ x + z, data = d)
  sml <- sqrt(sum(residuals(l)^2) / n)
  co <- summary(f)$coefficients
  db <- abs(unname(f$estimates$beta) - unname(coef(l))) /
    co$mu[, "Std. Error"]
  ds <- abs(unname(f$estimates$betad[1]) - log(sml)) /
    co$sigma[, "Std. Error"]
  llex <- sum(dnorm(d$y, fitted(l), sml, log = TRUE))
  X[[length(X) + 1L]] <- data.frame(
    build = lab, sd = sdv, seed = s, code = f$opt$convergence,
    beta_se = max(db), sigma_se = ds,
    llgap = llex - as.numeric(logLik(f)),
    unit = if (is.null(f$par_units)) 1 else max(f$par_units))
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev-log/lmprec-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
for (sdv in unique(X$sd)) {
  x <- X[X$sd == sdv, ]
  cat(sprintf(paste0("%s sd %-6g | median max|beta - lm|/SE %.2e (max %.2e) |",
                     " median |log sigma - ML|/SE %.2e | max ll gap %.2e |",
                     " code!=0 %d\n"),
              lab, sdv, median(x$beta_se), max(x$beta_se),
              median(x$sigma_se), max(x$llgap), sum(x$code != 0)))
}
