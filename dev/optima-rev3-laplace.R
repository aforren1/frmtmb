# Reviewer of lane optima, final check: check_laplace() after punch
# round 2, on designs the lane did not run: a MAP fit whose prior moves
# the estimate (so a dropped prior would show), two mo() terms, mo()
# with (1 | g) through laplace = TRUE, and mo() in sigma.
#   Rscript dev/optima-rev3-laplace.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
run <- function(label, fit, ...) {
  msg <- character()
  r <- withCallingHandlers(
    tryCatch(check_laplace(fit, chains = 2, iter = 2000, seed = 5,
                           refresh = 0, cores = 1, ...),
             error = function(e) e),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msg <<- c(msg, paste("WARN", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  cat("\n==", label, "==\n")
  for (m in unique(msg)) {
    cat("  msg:", substr(gsub("\n", " ", m), 1, 230), "\n")
  }
  if (inherits(r, "error")) {
    cat("  ERROR", conditionMessage(r), "\n")
  } else {
    print(format(r, digits = 3), row.names = FALSE)
  }
  invisible(r)
}
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
d$age <- rnorm(100, mean = 40, sd = 10)
# a prior that moves bsp_moincome from 15.4: if the retape dropped the
# fit's own prior, its post_mean would sit near 15.4, far from ml
fit_map <- frm(bf(ls ~ mo(income)), data = d, family = gaussian(),
               prior = set_prior("normal(0, 1)", class = "b"))
cat("MAP estimate of bsp:", fixef(fit_map)["moincome", "Estimate"], "\n")
run("MAP fit1, b ~ normal(0, 1)", fit_map)
run("two mo() terms, ls ~ mo(income) * age", frm(ls ~ mo(income) * age,
                                                 data = d))
set.seed(5)
n <- 300
dg <- data.frame(x1 = sample(0:3, n, TRUE), x2 = sample(0:4, n, TRUE),
                 z = rnorm(n), g = factor(sample(1:15, n, TRUE)))
dg$y <- c(0, 1, 1, 2)[dg$x1 + 1] + c(0, 0.3, 0.6, 0.6, 1)[dg$x2 + 1] +
  rnorm(15, sd = 0.5)[dg$g] + rnorm(n)
fg <- frm(bf(y ~ mo(x1) + mo(x2) + (1 | g)), data = dg, family = gaussian())
run("mo() + (1 | g), laplace = TRUE", fg, laplace = TRUE)
run("mo() + (1 | g), full", fg)
dg$ys <- 0.3 * dg$z + rnorm(n, sd = exp(c(0, 0.3, 0.3, 0.6)[dg$x1 + 1]))
run("mo() in sigma", frm(bf(ys ~ z, sigma ~ mo(x1)), data = dg,
                         family = gaussian()))
