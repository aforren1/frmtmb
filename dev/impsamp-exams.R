# Exams_2 against one-dimensional integrate(): is the corrected logLik right?
# Log dev/impsamp-out/exams.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
d <- data.frame(k = c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36,
                      39, 36, 35),
                n = 40L, zero = 0, id = factor(1:15))
fo <- bf(k | trials(n) ~ 1 + (1 | id), mu2 ~ 0 + offset(zero),
         theta1 ~ 0 + offset(zero))
imp <- suppressWarnings(suppressMessages(
  frm(fo, family = mixture(binomial, binomial), data = d, importance = 512)))
print(imp$opt$par)
b <- imp$opt$par[names(imp$opt$par) == "beta"]
s <- exp(imp$opt$par[names(imp$opt$par) == "theta"])
ex <- function(b, s) {
  sum(vapply(seq_len(nrow(d)), function(i) {
    f <- function(u) {
      (0.5 * dbinom(d$k[i], 40, plogis(b + u)) +
         0.5 * dbinom(d$k[i], 40, 0.5)) * dnorm(u, 0, s)
    }
    log(integrate(f, -8 * s, 8 * s, subdivisions = 500L,
                  rel.tol = 1e-10)$value)
  }, 0))
}
cat("corrected", format(as.numeric(logLik(imp)), digits = 9), "mcse",
    format(imp$importance$mcse, digits = 3), "exact at its estimates",
    format(ex(b, s), digits = 9), "\n")
lap <- suppressWarnings(suppressMessages(
  frm(fo, family = mixture(binomial, binomial), data = d)))
cat("Laplace", format(as.numeric(logLik(lap)), digits = 9), "exact at its",
    "estimates", format(ex(lap$opt$par[names(lap$opt$par) == "beta"],
                           exp(lap$opt$par[names(lap$opt$par) == "theta"])),
                        digits = 9), "\n")
