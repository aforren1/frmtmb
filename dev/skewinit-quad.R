# The escape runs on the quadrature tape too (it is skipped only under
# importance=). Checks that it does not break the quadrature path, and
# that a forced stall is still rescued there.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

set.seed(1)
n <- 200
g <- factor(rep(1:20, length.out = n))
xs <- -abs(rnorm(n)) * 3
y <- xs + rnorm(20, 0, 0.6)[as.integer(g)] +
  (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dd <- data.frame(y = y, xs = xs, g = g)
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
m3 <- skew(y)
bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))

show <- function(lab, ...) {
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd, ...))
  e <- f$opt[["stationary_escape"]]
  cat(sprintf("%-28s ll %16.9f alpha %10.4f esc %s\n", lab,
              as.numeric(logLik(f)),
              f$estimates$betad[["alpha_(Intercept)"]],
              if (is.null(e)) "-" else
                paste0(e[["starts"]], "/", format(e[["gain"]], digits = 5))))
  invisible(f)
}
cat("== Laplace\n")
show("default")
show("raw-skew start", start = bad)
cat("\n== importance = 200 (escape deliberately skipped)\n")
show("default", importance = 200)
cat("\n== quadrature = TRUE\n")
show("default", quadrature = TRUE)
show("raw-skew start", quadrature = TRUE, start = bad)
