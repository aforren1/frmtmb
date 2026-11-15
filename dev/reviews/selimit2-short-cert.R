# Short review of lane selimit2's fix round, F2: se_limit_joint_cert() on
# designs near its boundary. Prints each fit's verdict and SEs.
#   Rscript selimit2-short-cert.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
show <- function(label, expr) {
  ws <- character(0)
  ms <- character(0)
  fit <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(fit)))),
                 error = function(e) NA)
  li <- fit$cache$se_limit
  cat(sprintf("\n== %s | code %d\n   par: %s\n   se: %s\n", label,
              fit$opt$convergence,
              paste(format(fit$opt$par, digits = 5), collapse = " "),
              paste(format(se, digits = 4), collapse = " ")))
  if (!is.null(li)) {
    cat("   limit:", paste(names(li$kind), li$kind, collapse = "; "),
        "| joint", paste(li$joint, collapse = ","), "\n")
  }
  for (w in ws) cat("   W:", substr(w, 1, 200), "\n")
  for (m in ms) cat("   M:", substr(m, 1, 200), "\n")
}
# a logistic separated at x = 10 (shifted, so only the joint direction
# runs), with a row exactly at 10 (quasi-complete: MLE infinite), and with
# a row 1e-4 on the wrong side (still separable by moving the boundary:
# infinite), and a reversed pair 1e-4 apart (finite)
x0 <- c(5:9, 11:15)
y0 <- c(rep(0, 5), rep(1, 5))
mk <- function(x, y) data.frame(x = rep(x, 3), y = rep(y, 3))
show("bernoulli quasi, a row at the boundary",
     frm(y ~ x, data = mk(c(x0, 10, 10), c(y0, 0, 1)), family = bernoulli()))
show("bernoulli, a 1 just below the boundary",
     frm(y ~ x, data = mk(c(x0, 10 - 1e-4), c(y0, 1)), family = bernoulli()))
show("bernoulli, reversed pair 1e-4 apart (finite)",
     frm(y ~ x, data = mk(c(x0, 10 - 5e-5, 10 + 5e-5), c(y0, 1, 0)),
         family = bernoulli()))
# poisson: a level with one count of 1 among 200 zeros (finite), and the
# same with every count 0 (infinite)
set.seed(8)
d <- data.frame(g = factor(rep(c("a", "b"), c(200, 50))))
d$y <- c(replace(integer(200), 17, 1L), stats::rpois(50, 3))
show("poisson y ~ g, level a one count of 200", frm(y ~ g, data = d,
                                                     family = poisson()))
d$y[17] <- 0L
show("poisson y ~ g, level a all zero", frm(y ~ g, data = d,
                                             family = poisson()))
# poisson y ~ x, one count of 1 in the middle of the zeros, one at the top
x <- seq(0, 10, length.out = 40)
d2 <- data.frame(x = x, y = c(rep(0, 19), 1, rep(0, 19), 1))
show("poisson y ~ x, a count mid-range and one at the top (finite)",
     frm(y ~ x, data = d2, family = poisson()))
# zi ~ f, level a with 3 nonzero rows of 40, and zi ~ x with the nonzero
# rows only above x = 0 but one at x = -2.9
set.seed(2)
d3 <- data.frame(f = factor(rep(c("a", "b"), each = 40)))
d3$y <- c(c(rep(0, 37), 2, 1, 4), stats::rpois(40, 3))
show("zip zi ~ f, level a 3 nonzero of 40",
     frm(bf(y ~ 1, zi ~ f), data = d3, family = zero_inflated_poisson()))
x4 <- seq(-3, 3, length.out = 80)
y4 <- ifelse(x4 < 0, 0L, stats::rpois(80, 3))
show("zip zi ~ x, every row below 0 a zero (zi to 1 there)",
     frm(bf(y ~ 1, zi ~ x), data = data.frame(x = x4, y = y4),
         family = zero_inflated_poisson()))
y4[1] <- 2L
show("zip zi ~ x, the same with a 2 at x = -3",
     frm(bf(y ~ 1, zi ~ x), data = data.frame(x = x4, y = y4),
         family = zero_inflated_poisson()))
