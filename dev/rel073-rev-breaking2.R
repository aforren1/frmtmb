# Release reviewer, 0.73.0: two follow-ups of dev/rel073-rev-breaking.R:
# a valid by-split block under quadrature (the NEWS says it is unchanged),
# and a fractional trunc() bound on a count (the NEWS says it is refused).
#   Rscript dev/rel073-rev-breaking2.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib: frmtmb", format(packageVersion("frmtmb")), "\n")
set.seed(5)
nq <- 240
dq <- data.frame(x = rnorm(nq), ga = factor(rep(1:20, each = 12)))
dq$f <- factor(ifelse(as.integer(dq$ga) <= 10, "u", "v"))
dq$y <- rbinom(nq, 1, plogis(0.3 * dq$x + rnorm(20, 0, 0.8)[dq$ga]))
r <- tryCatch(suppressWarnings(frm(y ~ x + (1 | gr(ga, by = f)),
                                   family = bernoulli(), data = dq,
                                   quadrature = TRUE)),
              error = function(e) conditionMessage(e))
cat("by-split quadrature:", if (is.character(r)) paste("ERROR", substr(r, 1, 200))
    else sprintf("logLik %.6f", as.numeric(logLik(r))), "\n")
d <- data.frame(x = rnorm(100)); d$k <- rpois(100, 4)
for (lbv in c(1.5, 2)) {
  dd <- d[d$k >= 2, ]; dd$lb <- lbv
  r <- tryCatch(suppressWarnings(frm(k | trunc(lb = lb) ~ x, family = poisson(),
                                     data = dd)),
                error = function(e) conditionMessage(e))
  cat("trunc(lb =", lbv, "):", if (is.character(r)) paste("ERROR", substr(r, 1, 200))
      else sprintf("logLik %.6f", as.numeric(logLik(r))), "\n")
}
